import torch
import torch.nn as nn
from torch.nn.utils import clip_grad_norm_
from torch.utils.data import Dataset, DataLoader, random_split
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
from transformers import AutoTokenizer, AutoModel
import pandas as pd
from PRP49_ban import BANLayer
from torch.nn.utils.parametrizations import weight_norm
import copy
import numpy as np

# -------------------- 模型定义 --------------------
class BindingEnergyPredictor(nn.Module):
    def __init__(self,
                 esm_model_name='facebook/esm2_t6_8M_UR50D',
                 lassoesm_model_name='./LassoESM',
                 h_dim=1024,          # BAN 隐藏维度
                 n_heads=8,           # 注意力头数
                 k=3,                 # 低秩分解参数
                 dropout=0.3,
                 ban_dropout=0.3,
                 freeze_esm=False,    # ESM 冻结
                 lassoesm_unfreeze=0, # LassoESM 解冻，传入 -1 解冻全部参数
                 output_dim=1):
        super().__init__()

        # ESM编码器
        self.esm = AutoModel.from_pretrained(esm_model_name)
        self.esm_emb_dim = self.esm.config.hidden_size  # 320

        # LassoESM编码器
        self.lassoesm = AutoModel.from_pretrained(lassoesm_model_name)
        self.lassoesm_emb_dim = self.lassoesm.config.hidden_size  # 1280

        if freeze_esm:
            for param in self.esm.parameters():
                param.requires_grad = False

        for param in self.lassoesm.parameters():
            param.requires_grad = False

        if lassoesm_unfreeze == -1:
            for param in self.lassoesm.parameters():
                param.requires_grad = True
        elif lassoesm_unfreeze > 0:
            total_layers = len(self.lassoesm.encoder.layer)
            for i in range(total_layers - lassoesm_unfreeze, total_layers):
                for param in self.lassoesm.encoder.layer[i].parameters():
                    param.requires_grad = True

        # 双线性注意力融合层（v=protein特征， q=lasso特征）
        self.ban = weight_norm(
            BANLayer(v_dim=self.esm_emb_dim,
                     q_dim=self.lassoesm_emb_dim,
                     h_dim=h_dim,
                     h_out=n_heads,
                     k=k,
                     dropout=ban_dropout,),
            name='h_mat', dim=None
        )

        # 回归头
        self.fc1 = nn.Linear(h_dim, 512)
        self.fc2 = nn.Linear(512, 128)
        self.out = nn.Linear(128, output_dim)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

    def extract_aa_tokens(self, token_embeddings, input_ids):
        """
        去掉 <CLS>/<EOS>/<PAD>，只保留氨基酸 token（id >= 4）
        """
        aa_mask = input_ids >= 4
        aa_embs = []
        for i in range(token_embeddings.size(0)):
            aa_embs.append(token_embeddings[i][aa_mask[i]])
        max_len = max(emb.size(0) for emb in aa_embs) if aa_embs else 1
        dim = token_embeddings.size(-1)   # 自动适应编码器维度
        padded = torch.zeros(len(aa_embs), max_len, dim,
                             device=token_embeddings.device)
        for i, emb in enumerate(aa_embs):
            padded[i, :emb.size(0)] = emb
        return padded

    def forward(self, seq1_ids, seq1_mask, seq2_ids, seq2_mask):
        out1 = self.esm(input_ids=seq1_ids, attention_mask=seq1_mask)
        out2 = self.lassoesm(input_ids=seq2_ids, attention_mask=seq2_mask)

        feat1 = self.extract_aa_tokens(out1.last_hidden_state, seq1_ids)   # (B, L1, 320)
        feat2 = self.extract_aa_tokens(out2.last_hidden_state, seq2_ids)   # (B, L2, 1280)

        # 双线性注意力融合
        fused, _ = self.ban(feat1, feat2)   # (B, h_dim)

        # MLP 预测
        x = self.fc1(fused)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.out(x)
        return x.squeeze(-1)


# -------------------- 数据集 --------------------
class BindingEnergyDataset(Dataset):
    def __init__(self, csv_path, max_len=512):
        self.df = pd.read_csv(csv_path)
        self.max_len = max_len

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        protein_seq = row['protein seq']
        lasso_seq = row['lasso seq']
        energy = row['energy']
        return {
            'protein_seq': protein_seq,
            'lasso_seq': lasso_seq,
            'energy': torch.tensor(energy, dtype=torch.float32)
        }


# -------------------- 批处理函数 --------------------
def collate_fn(batch, tokenizer, max_len):
    protein_seqs = [item['protein_seq'] for item in batch]
    lasso_seqs = [item['lasso_seq'] for item in batch]
    energies = torch.stack([item['energy'] for item in batch], dim=0)

    tok1 = tokenizer(protein_seqs, return_tensors='pt',
                     max_length=max_len, truncation=True, padding=True)
    tok2 = tokenizer(lasso_seqs, return_tensors='pt',
                     max_length=max_len, truncation=True, padding=True)

    return {
        'seq1_ids': tok1['input_ids'],
        'seq1_mask': tok1['attention_mask'],
        'seq2_ids': tok2['input_ids'],
        'seq2_mask': tok2['attention_mask'],
        'energy': energies
    }

# -------------------- 随机种子 --------------------
def set_seed(seed=42):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# -------------------- 训练主程序 --------------------
def main():
    set_seed(42)
    csv_path = r"C:\autodockvina\prp49\outdata\out.csv"
    esm_model_name = "facebook/esm2_t6_8M_UR50D"
    lassoesm_path = "./LassoESM"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    batch_size = 32
    epochs = 500
    lr = 8e-5
    max_len = 1024
    train_ratio = 0.8
    patience = 150
    min_lr = 0.5e-6
    dropout = 0.5
    ban_dropout = 0.4
    T_0 = 20
    T_mult = 2

    # ========== 模型、优化算法、损失函数、学习率调度器 ==========
    model = BindingEnergyPredictor(
        esm_model_name=esm_model_name,
        lassoesm_model_name=lassoesm_path,
        h_dim=256,
        n_heads=3,
        k=2,
        dropout=dropout,
        ban_dropout=ban_dropout,
        freeze_esm=False,
        lassoesm_unfreeze=1
    )
    model.to(device)
    model.train()

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    loss_fn = nn.MSELoss()

    scheduler = CosineAnnealingWarmRestarts(
        optimizer,
        T_0=T_0,
        T_mult=T_mult,
        eta_min=min_lr,
    )

    # ========== 加载 tokenizer 和数据集 ==========
    tokenizer = AutoTokenizer.from_pretrained(esm_model_name)
    full_dataset = BindingEnergyDataset(csv_path, max_len=max_len)

    train_size = int(train_ratio * len(full_dataset))
    print("训练集大小：", train_size)
    val_size = len(full_dataset) - train_size
    train_dataset, val_dataset = random_split(
        full_dataset,
        [train_size, val_size],
        generator=torch.Generator().manual_seed(42)
    )

    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True,
        collate_fn=lambda batch: collate_fn(batch, tokenizer, max_len)
    )
    val_loader = DataLoader(
        val_dataset, batch_size=batch_size, shuffle=False,
        collate_fn=lambda batch: collate_fn(batch, tokenizer, max_len)
    )


    # ========== 早停机制 ==========
#    best_val_loss = float('inf')
    best_r2 = -float('inf')
    best_epoch = 0
    early_stop_counter = 0
    best_model_state = None

    # ========== 训练循环 ==========
    print("start training")
    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        for batch in train_loader:
            seq1_ids = batch['seq1_ids'].to(device)
            seq1_mask = batch['seq1_mask'].to(device)
            seq2_ids = batch['seq2_ids'].to(device)
            seq2_mask = batch['seq2_mask'].to(device)
            energy = batch['energy'].to(device)

            optimizer.zero_grad()
            pred = model(seq1_ids, seq1_mask, seq2_ids, seq2_mask)
            loss = loss_fn(pred, energy)
            loss.backward()
            clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            train_loss += loss.item()

        avg_train_loss = train_loss / len(train_loader)

        model.eval()
        val_loss = 0.0
        all_preds = []
        all_labels = []
        with torch.no_grad():
            for batch in val_loader:
                seq1_ids = batch['seq1_ids'].to(device)
                seq1_mask = batch['seq1_mask'].to(device)
                seq2_ids = batch['seq2_ids'].to(device)
                seq2_mask = batch['seq2_mask'].to(device)
                energy = batch['energy'].to(device)

                pred = model(seq1_ids, seq1_mask, seq2_ids, seq2_mask)
                loss = loss_fn(pred, energy)
                val_loss += loss.item()
                all_preds.extend(pred.cpu().numpy().tolist())
                all_labels.extend(energy.cpu().numpy().tolist())

        avg_val_loss = val_loss / len(val_loader)

        all_preds = np.array(all_preds)
        all_labels = np.array(all_labels)
        rmse = np.sqrt(np.mean((all_preds - all_labels) ** 2))
        mae = np.mean(np.abs(all_preds - all_labels))
        ss_res = np.sum((all_labels - all_preds) ** 2)
        ss_tot = np.sum((all_labels - np.mean(all_labels)) ** 2)
        ss_tot = float(ss_tot)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        scheduler.step()

        current_lr = optimizer.param_groups[0]['lr']

        print(f"Epoch {epoch + 1:3d} | Train Loss: {avg_train_loss:.6f} | "
              f"Val Loss: {avg_val_loss:.6f} | RMSE: {rmse:.4f} | MAE: {mae:.4f} | R2: {r2:.4f} | LR: {current_lr:.2e}")

#        if avg_val_loss < best_val_loss:
#            best_val_loss = avg_val_loss
#            best_epoch = epoch + 1
#            early_stop_counter = 0
#            best_model_state = copy.deepcopy(model.state_dict())
#        else:
#            early_stop_counter += 1
#            if early_stop_counter >= patience:
#                print(f"早停触发！验证 loss 连续 {patience} 轮未改善，训练在第 {epoch+1} 轮停止。")
#                break

        if r2 > best_r2 - 0.05:
            if r2 > best_r2:
                best_r2 = r2
                best_epoch = epoch + 1
                best_model_state = copy.deepcopy(model.state_dict())
            early_stop_counter = 0
        else:
            early_stop_counter += 1
            if early_stop_counter >= patience:
                print("早停触发！")
                break

    #            best_epoch = epoch + 1
    #           early_stop_counter = 0
    #            best_model_state = copy.deepcopy(model.state_dict())
    #        else:
    #            early_stop_counter += 1
    #            if early_stop_counter >= patience:
    #                print(f"早停触发！验证 loss 连续 {patience} 轮未改善，训练在第 {epoch+1} 轮停止。")
    #                break

    if best_model_state is not None:
        model.load_state_dict(best_model_state)
        print(f"已加载最佳模型（Epoch {best_epoch}，Best R2: {best_r2:.6f}）")
        torch.save(best_model_state, "lpiLasso.pth")
    else:
        torch.save(model.state_dict(), "lpiLasso.pth")

    print("training finished")

if __name__ == "__main__":
    main()