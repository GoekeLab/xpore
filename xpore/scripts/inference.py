import argparse
import torch
import pandas as pd
import numpy as np
import os
import h5py
from torch.utils.data import DataLoader, Dataset
from itertools import product


class MultiInstanceNNEmbedding(torch.nn.Module):
    def __init__(self, dim_cov, p=1, embedding_dim=2):
        """
        fix_classifier: if True the last dense layer is fixed to a random orthogonal matrix
        """
        super(MultiInstanceNNEmbedding, self).__init__()
        self.p = p
        self.embedding_dim = embedding_dim
        self.embedding = torch.nn.Embedding(18, self.embedding_dim)
        self.linear1 = torch.nn.Linear(dim_cov + self.embedding_dim, 150, bias=True)
        self.linear2 = torch.nn.Linear(150, p, bias=True)
        self.linear3 = torch.nn.Linear(5 * p, 150, bias=True)
        self.linear4 = torch.nn.Linear(150, 1, bias=True)

    def forward(self, x, kmer, indices):
        """ compute probability at site level"""
        kmer_embedding = self.embedding(kmer)
        x = torch.cat([x, kmer_embedding], axis=1)
        x = torch.relu(self.linear1(x))
        x = torch.relu(self.linear2(x))
        x = self.aggregate(x, indices)
        x = x.view(-1, 5 * self.p)
        x = torch.relu(self.linear3(x))
        out = torch.sigmoid(self.linear4(x))
        return out

    def aggregate(self, x, indices):
        grouped_tensors = [x[idx] for idx in indices]
        mean = torch.stack([torch.mean(tensor, axis=0) for tensor in grouped_tensors])
        std = torch.stack([torch.var(tensor, axis=0) for tensor in grouped_tensors])
        maximum = torch.stack([torch.max(tensor, axis=0).values for tensor in grouped_tensors])
        minimum = torch.stack([torch.min(tensor, axis=0).values for tensor in grouped_tensors])
        median = torch.stack([torch.median(tensor, axis=0).values for tensor in grouped_tensors])
        aggregated_tensors = torch.stack([mean, std, minimum, median, maximum], axis=1)
        return aggregated_tensors


class ValDS(Dataset):

    def __init__(self, norm_constant, data_dir, sites=None):
        self.data_dir = data_dir
        self.norm_constant = norm_constant.set_index("0")
        
        self.all_kmers = list(["".join(x) for x in product(['A', 'G', 'T'], ['G', 'A'], ['A'], ['C'], ['A', 'C', 'T'])])
        self.kmer_to_int = {self.all_kmers[i]: i for i in range(len(self.all_kmers))}
        self.int_to_kmer =  {i: self.all_kmers[i] for i in range(len(self.all_kmers))}
        
        if sites is None:
            all_files = os.listdir(data_dir)
            self.sites = np.array([os.path.join(data_dir, fname) for fname in all_files])
            self.kmers = np.array([self.kmer_to_int[x.split("_")[2]] for x in all_files])

        else:
            self.sites = np.array([os.path.join(data_dir, fname) for fname in sites])
            self.kmers = np.array([self.kmer_to_int[x.split("_")[2]] for x in sites])
            
            
    def __len__(self):
        return len(self.sites)

    def __getitem__(self, idx):
        kmer = self.kmers[idx]
        norm_info = self.norm_constant.loc[self.int_to_kmer[kmer]].values
        mean, std = norm_info[:3], norm_info[3:]
        if torch.is_tensor(idx):
            idx = idx.tolist()
        f = h5py.File(self.sites[idx], 'r')
        X = (f['X'][:] - mean) / std
        n_reads = len(X)
        f.close()
        return (torch.Tensor(X),
                torch.LongTensor([kmer]).repeat(len(X)),
                n_reads)

                
def custom_collate_val(batch):
    return (torch.cat([item[0] for item in batch]),
            torch.cat([item[1] for item in batch]),
            assign_group_to_index_val(batch))


def assign_group_to_index_val(batch):
    curr_idx = 0
    idx_per_group = []
    for i in range(len(batch)):
        num_reads = batch[i][2]
        idx_per_group.append(np.arange(curr_idx, curr_idx + num_reads))
        curr_idx += num_reads
    return np.array(idx_per_group)


def extract_positions(tx_dir):
    fnames = os.listdir(tx_dir)
    fpaths = [os.path.join(tx_dir, fname) for fname in fnames]
    transcripts = [x.split("_")[0] for x in fnames]
    positions = [np.array(x.split("_"))[1] for x in fnames]
    n_samples = [int(x.split("_")[-1].split(".hdf5")[0]) for x in fnames]
    kmers = [x.split("_")[2] for x in fnames]
    return pd.DataFrame({'filepath': fpaths, 'position': positions, 'transcript_id': transcripts,
                         'n_samples': n_samples,
                         'kmer': kmers})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="a script to preprocess all files in a nanopolish event align directory")
    parser.add_argument('-i', '--input', dest='input_dir', default=None,
                        help="Directory containing the inference folder to predict on")
    parser.add_argument('-o', '--output', dest='out_dir', default=None,
                        help="Save directory for the prediction results")
    parser.add_argument('-m', '--model', dest='model_path', default=None,
                        help="Path to directory containing norm constant and state dictionary for torch model")
    parser.add_argument('-d', '--device', dest='device', default='cpu',
                        help="cpu or cuda to run the inference on")
    parser.add_argument('-n', '--n_processors', dest='n_processors', default=None,
                        help="number of workers for dataloader")
    args = parser.parse_args()
    state_dict = torch.load(os.path.join(args.model_path, "best_model.pt"))
    norm_constant = pd.read_csv(os.path.join(args.model_path, "norm_constant.csv"))    
    device = args.device
    model = MultiInstanceNNEmbedding(dim_cov=3, p=8, embedding_dim=2).to(device)
    model.load_state_dict(state_dict)
    
    data_dir = os.path.join(args.input_dir, "inference")
    df = extract_positions(data_dir)
    ds = ValDS(norm_constant, data_dir)
    dl = DataLoader(ds, num_workers=int(args.n_processors), batch_size=150,
                    shuffle=False, collate_fn=custom_collate_val)
    y_preds = []
    
    for _, inp in enumerate(dl):
        out = model(inp[0].to(device), inp[1].to(device), inp[2])
        y_preds.extend(out.detach().cpu().numpy())

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    df["proba"] = np.array(y_preds)
    df[["transcript_id", "kmer", "n_samples", "proba", "filepath"]]\
        .to_csv(os.path.join(args.out_dir, "m6Anet_predictions.csv.gz"), index=False)
