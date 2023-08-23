import argparse
import hashlib
import re
import os
import numpy as np
from jax.tree_util import tree_map
from scipy.special import softmax
from string import ascii_uppercase, ascii_lowercase
from Bio import SeqIO
import torch
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Run ESMFold")
    parser.add_argument("fasta_file", help="Path to input fasta file")
    args = parser.parse_args()
    return args


def parse_output(output):
    pae = (output["aligned_confidence_probs"][0] * np.arange(64)).mean(-1) * 31
    plddt = output["plddt"][0, :, 1] * 100
    # printing parsing output data
    #print("Parsed Output: ", pae, plddt)

    bins = np.append(0, np.linspace(2.3125, 21.6875, 63))
    sm_contacts = softmax(output["distogram_logits"], -1)[0]
    sm_contacts = sm_contacts[..., bins < 8].sum(-1)
    xyz = output["positions"][-1, 0, :, 1]
    mask = output["atom37_atom_exists"][0, :, 1] == 1
    o = {
        "pae": pae[mask, :][:, mask],
        "plddt": plddt[mask],
        "sm_contacts": sm_contacts[mask, :][:, mask],
        "xyz": xyz[mask],
    }
    # if "contacts" in output["lm_output"]:
    #    lm_contacts = output["lm_output"]["contacts"].astype(float)[0]
    #    o["lm_contacts"] = lm_contacts[mask, :][:, mask]
    return o


def get_hash(x):
    return hashlib.sha1(x.encode()).hexdigest()


def main():
    args = parse_args()
    alphabet_list = list(ascii_uppercase + ascii_lowercase)
    num_recycles = 3
    get_LM_contacts = False
    copies = 1
    chain_linker = 25
    samples = None
    masking_rate = 0.15
    stochastic_mode = "LM"
    thaliania = True

    print("Initial arguments: ", args)

    if "model" not in dir():
        print("Loading model...")
        model = torch.load("esmfold.model")
        model.cuda().requires_grad_(False)

    for record in tqdm(SeqIO.parse(args.fasta_file, "fasta")):
        jobname = record.id
        sequence = str(record.seq)
        # sequence preprocessing
        sequence = re.sub("[^A-Z:]", "", sequence.replace("/", ":").upper())
        sequence = re.sub(":+", ":", sequence)
        sequence = re.sub("^[:]+", "", sequence)
        sequence = re.sub("[:]+$", "", sequence)
        sequence = ":".join([sequence] * copies)

        if thaliania == True:
            ID = jobname
        else:
            ID = jobname + "_" + get_hash(sequence)[:5]
        seqs = sequence.split(":")
        lengths = [len(s) for s in seqs]
        length = sum(lengths)
        #print("length", length)

        u_seqs = list(set(seqs))
        if len(seqs) == 1:
            mode = "mono"
        elif len(u_seqs) == 1:
            mode = "homo"
        else:
            mode = "hetero"

        #print("Job Mode: ", mode)

        if length > 700:
            model.trunk.set_chunk_size(64)
        else:
            model.trunk.set_chunk_size(128)

        best_pdb_str = None
        best_ptm = 0
        best_output = None
        traj = []

        num_samples = 1 if samples is None else samples
        for seed in range(num_samples):
            torch.cuda.empty_cache()
            if samples is None:
                seed = "default"
                mask_rate = 0.0
                model.train(False)
            else:
                torch.manual_seed(seed)
                mask_rate = masking_rate if "LM" in stochastic_mode else 0.0
                model.train("SM" in stochastic_mode)

            #print("Running model inference...")
            output = model.infer(
                sequence,
                num_recycles=num_recycles + 1,
                chain_linker="X" * chain_linker,
                residue_index_offset=512,
                # mask_rate=mask_rate,
                # return_contacts=get_LM_contacts,
            )

            pdb_str = model.output_to_pdb(output)[0]
            output = tree_map(lambda x: x.cpu().numpy(), output)
            ptm = output["ptm"]
            plddt = output["plddt"][..., 1].mean()
            traj.append(parse_output(output))
            # print(f"{seed} ptm: {ptm:.3f} plddt: {plddt:.3f}")
            if ptm > best_ptm:
                best_pdb_str = pdb_str
                best_ptm = ptm
                best_output = output
            # os.system(f"mkdir -p {ID}")
            if samples is None:
                if thaliania == True:
                    pdb_filename = f"./proteins_data/{ID}/{ID}_ESM.pdb"
                else:
                    pdb_filename = f"{ID}/{ID}_ESM.pdb"
            else:
                pdb_filename = f"{ID}/ptm{ptm:.3f}_r{num_recycles}_seed{seed}_{stochastic_mode}_mr{masking_rate:.2f}.pdb"

            with open(pdb_filename, "w") as out:
                out.write(pdb_str)
                #print("Output written to: ", pdb_filename)


if __name__ == "__main__":
    main()
