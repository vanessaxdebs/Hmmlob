import os
import subprocess
import random
import sys, yaml, requests

from Bio import AlignIO
from Bio.PDB import PDBList, PDBParser, PDBIO, Structure, Model, Chain

import seaborn as sns


UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
PDBE_MAPPINGS_BASE_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/interpro"
EBI_FASTA_BASE_URL = "https://www.ebi.ac.uk/pdbe/entry/pdb"
RCSB_BASE_URL_FOR_ALL_PDBS="https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
UNIPROT_HEADERS = {"accept": "application/json"}
OUTPUT_LINE_TEMPLATE_ALL_STRUCTURES="Initial search space (Z)"
OUTPUT_LINE_TEMPLATE_PASSED_STRUCTURES="Domain search space  (domZ)"

def build_query(domain, limit):
    """Build a query for UniProt API.
    """
    return {
        "query": f"(xref:interpro-{domain}) AND reviewed:true AND (database:pdb)",
        "fields": [
            "accession",
            "xref_pdb"
        ],
        "sort": "accession desc",
        "size": f"{limit}"
    }


def do_request(url, headers=None, params=None):
    """
    Do HTTP GET request.
    :param url: URL
    :param headers: Dictionary of HTTP Headers
    :param params: Dictionary of query parameters
    :return: Response object
    """
    response = requests.get(
        url,
        headers=headers,
        params=params,
    )

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    return response


def get_ids(domain, limit):
    """
    Obtain PDB IDs that contain the domain.
    :param domain: the domain for a search query
    :param limit: the limit for a search query
    :return: the list of PDB IDs
    """
    data = do_request(
        UNIPROT_BASE_URL,
        UNIPROT_HEADERS,
        build_query(domain, limit),
    ).json()

    pdb_ids = []
    for uniprot in data["results"]:
        if (len(uniprot["uniProtKBCrossReferences"]) != 0 and
                uniprot["uniProtKBCrossReferences"][0]["database"] == "PDB"):
            pdb_ids.append(uniprot["uniProtKBCrossReferences"][0]["id"])

    return pdb_ids


def get_chains(domain, pdb_id):
    """
    Obtain a set of chains that contain the domain.
    :param domain: the domain for a search query
    :param pdb_id: the PDB ID for a search query
    :return: the set of chains of the PDB file that contain this domain
    """
    try:
       resp = do_request(f"{PDBE_MAPPINGS_BASE_URL}/{pdb_id.lower()}")
       chains = []
       for found_domain, info in resp.json().get(pdb_id.lower(), {}).get("InterPro", {}).items():
           if found_domain != domain:
               continue
           for mapping in info.get("mappings", []):
               chains.append(mapping["chain_id"])
       return set(chains)
    except Exception as e:
        return set()


def filter_by_resolution(pdb_file, resolution_threshold):
    """
    Filter PDB file by its resolution.
    :param pdb_file: the PDB file to filter
    :param resolution_threshold: the threshold for the resolution
    :return: Boolean flag that says if PDB should be discarded
    """
    resolution_line = subprocess.run(
        f"grep \"^REMARK   2 RESOLUTION.\" {pdb_file}",
        capture_output=True,
        shell=True,
    )

    line = resolution_line.stdout.strip()
    if line:
        parts = line.split()
        if len(parts) > 3:
            try:
                resolution = float(parts[3])
                if resolution > resolution_threshold:
                    return True
            except ValueError:
                return False
    else:
        return False


def download_pdb_files(
        pdbl,
        parser,
        io,
        pdb_ids,
        domain,
        output_dir,
        pdb_filter,
):
    """
    Download chains in the PDB format and filter them according to
    constraints set in the config.
    :param pdbl:
    :param parser:
    :param io:
    :param pdb_ids: the IDs of the PDB files to download
    :param domain: the domain of the interest
    :param output_dir: the output directory
    :param pdb_filter: the dictionary with constraints
    :return:
    """
    os.makedirs(output_dir, exist_ok=True)
    for pdb_id in sorted(pdb_ids):
        chains = get_chains(domain, pdb_id)
        if len(chains) == 0:
            continue

        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=f"{output_dir}/pdb", file_format="pdb")

        # Skip all NMR structures
        is_nmr = False
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("EXPDTA") and "NMR" in line.upper():
                    is_nmr = True
                    break

        if is_nmr:
            continue

        # Filter the PDB file by resolution
        if filter_by_resolution(pdb_file, pdb_filter["resolution_threshold"]):
            continue

        struct = parser.get_structure(pdb_id, pdb_file)

        # Get the first model only
        model = next(struct.get_models())

        for chain_id in chains:
            if chain_id not in model:
                continue

            out_name = os.path.join(f"{output_dir}/pdb", f"{pdb_id}_{chain_id}.pdb")

            # Filter the chain by its length
            chain = model[chain_id]
            length = sum(1 for _ in chain)
            if length > pdb_filter["max_residues_count"]:
                continue

            new_structure = Structure.Structure(pdb_id)
            new_model = Model.Model(0)
            new_chain = Chain.Chain(chain_id)

            for residue in chain:
                new_chain.add(residue.copy())  # ensure clean copy

            new_model.add(new_chain)
            new_structure.add(new_model)

            io.set_structure(new_structure)
            io.save(out_name)


def align(output_dir):
    """
    Align the structures using MUSTANG tool (Multiple Structural Alignment).
    :param output_dir: the output directory
    :return:
    """
    os.makedirs(f"{output_dir}/alignment", exist_ok=True)

    with open(f"{output_dir}/alignment/list.txt", 'w') as list:
        list.write("> .\n")
        for file in os.listdir(f"{output_dir}/pdb"):
            if file.endswith(".pdb"):
                list.write(f"+{output_dir}/pdb/{file}\n")

    subprocess.run(f"mustang -f {output_dir}/alignment/list.txt -o {output_dir}/alignment/alignment -F fasta",
                   shell=True)


def fetch_fasta_by_pdb_id(pdb_id):
    """
    Fetch FASTA file by PDB ID.
    :param pdb_id: ID of a PDB file
    :return:
    """
    try:
        return do_request(f"{EBI_FASTA_BASE_URL}/{pdb_id.lower()}/fasta").text
    except Exception as e:
        raise Exception(f"Failed to fetch FASTA for {pdb_id}")


def download_validation_set(output_dir, ids, fraction, domain):
    """
    Download a set of PDB structures that are going to be used for the validation of
    the model. For the validation we use the same structures that we used to extract chains,
    but only the subset of them (set validation_set_fraction in the config).
    :param output_dir: the output directory
    :param ids: ids of PDB files to download
    :param fraction: the percentage of the PDB files to download
    :param domain: the domain of the interest
    :return: the amount of FASTA files downloaded
    """
    number = int(len(ids) * fraction)
    os.makedirs(f"{output_dir}/validation", exist_ok=True)
    download_fasta(
        f"{output_dir}/validation/{domain}.fa",
        random.sample(ids, number)
    )
    return number

def download_validation_random_set(output_dir, number):
    """
    Download a set of PDB structures that are going to be used for the validation of
    the model. For the validation we use the same amount of random structures as
    the amount of the structures with the domain.
    :param output_dir: the output directory
    :param number: number of PDB files to download
    :return:
    """
    download_fasta(
        f"{output_dir}/validation/random.fa",
        random.sample(do_request(RCSB_BASE_URL_FOR_ALL_PDBS).json(), number)
    )

def download_fasta(output_file, ids):
    """
    Download FASTA files by PDB IDs.
    :param ids: List of PDB IDs
    :param output_file: output file
    :return:
    """
    with open(output_file, 'w') as file:
        for id in ids:
            file.write(fetch_fasta_by_pdb_id(id))


def fasta_to_stockholm(output_dir):
    """
    Converts FASTA file to the Stockholm alignment file
    :param output_dir: the output directory
    :return:
    """
    # Read the FASTA alignment
    alignment = AlignIO.read(f"{output_dir}/alignment/alignment.afasta", "fasta")

    # Save in Stockholm format
    AlignIO.write(alignment, f"{output_dir}/alignment/alignment.sto", "stockholm")


def build_model(output_dir, domain):
    """
    Build a profile HMM using hmmbuild routine of the HHMER.
    :param output_dir: the output directory
    :param domain: the domain of the interest
    :return:
    """
    # convert to Stockholm format
    fasta_to_stockholm(output_dir)

    os.makedirs(f"{output_dir}/model", exist_ok=True)

    # build a model
    subprocess.run(
        f"hmmbuild {output_dir}/model/{domain}.hmm {output_dir}/alignment/alignment.sto",
        shell=True,
    )


def validate(output_dir, domain, output_name, e_value_threshold):
    """
    Validate the model using hmmsearch routine of the HHMER.
    :param output_dir: the output directory
    :param domain: the domain of the interest
    :param output_name: the name of the output file
    :return:
    """
    subprocess.run(
        f"hmmsearch -E {e_value_threshold} {output_dir}/model/{domain}.hmm {output_dir}/validation/{output_name}.fa > {output_dir}/validation/{output_name}.output",
        shell=True,
    )

def get_metrics(
        output_dir,
        domain,
        random_output_name,
):
    """
    Calculate metrics (True Positive, False Positive, True Negative, False Negative).
    :param output_dir: the output directory
    :param domain: the domain of the interest
    :param random_output_name: the name of the hmmsearch output file (with random proteins)
    :return:
    """
    TP, FP, TN, FN = 0,0,0,0
    try:
        total_count = extract_the_metric(
            output_dir,
            domain,
            OUTPUT_LINE_TEMPLATE_ALL_STRUCTURES,
        )
        TP = extract_the_metric(
            output_dir,
            domain,
            OUTPUT_LINE_TEMPLATE_PASSED_STRUCTURES,
        )
        FN = total_count - TP

        total_count_of_random = extract_the_metric(
            output_dir,
            random_output_name,
            OUTPUT_LINE_TEMPLATE_ALL_STRUCTURES,
        )
        FP = extract_the_metric(
            output_dir,
            random_output_name,
            OUTPUT_LINE_TEMPLATE_PASSED_STRUCTURES,
        )
        TN = total_count_of_random - FP
    except Exception as e:
        raise Exception("Failed to extract all metrics")

    return TP, FP, TN, FN

def extract_the_metric(
    output_dir,
    filename,
    line_beginning
):
    """
    Extract metrics from the file
    :param output_dir: the output directory
    :param filename: the name of the file
    :param line_beginning: the mask for identify the line
    :return:
    """
    grep_line = subprocess.run(
        f"grep \"^{line_beginning}:\" {output_dir}/validation/{filename}.output",
        capture_output=True,
        shell=True,
    )
    line = grep_line.stdout.strip()
    if grep_line.stdout.strip():
        parts = line.split()
        if len(parts) < 3:
            raise Exception("Corrupted output of hmmsearch")
        return int(parts[4])
    raise Exception("Corrupted output of hmmsearch")

def search(
    output_dir,
    e_value_threshold,
    domain,
    db_fraction,
    db=None,
):
    """
    Search for sufficient proteins in the DB using pHMM.
    :param output_dir: the output directory
    :param e_value_threshold: E-value threshold
    :param domain: the domain of interest
    :param db_fraction: percentage of DB to use
    :param db: pass to the DB file
    :return:
    """
    os.makedirs(f"{output_dir}/search", exist_ok=True)
    # get all PDB IDs
    resp = do_request(RCSB_BASE_URL_FOR_ALL_PDBS).json()
    if db is None or db == "":
        download_fasta(
            f"{output_dir}/search/{domain}_candidates.fa",
            random.sample(resp, int(float(db_fraction)*len(resp)))
        )
        db = f"{output_dir}/search/{domain}_candidates.fa"
    subprocess.run(
        f"hmmsearch -E {e_value_threshold} {output_dir}/model/{domain}.hmm {db} > {output_dir}/search/{domain}.output",
        shell=True,
    )
    subprocess.run(
        "awk '$1 ~ /^[0-9]/ && $1+0 < 10 { print $9 }'"+ f" {output_dir}/search/{domain}.output | grep 'pdb' > {output_dir}/search/{domain}_ids",
        shell=True,
    )


if __name__ == "__main__":
    # Init biopython objects
    pdbl = PDBList()
    parser = PDBParser()
    io = PDBIO()

    # Read config
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    random.seed(config['seed'])

    # Get PDB identifiers from Uniprot
    pdb_ids = get_ids(
        config['target']['interpro_id'],
        config['target']['structures_number_limit'],
    )

    # Download PDB files for chains that contain the domain
    download_pdb_files(
        pdbl,
        parser,
        io,
        pdb_ids,
        config['target']['interpro_id'],
        config['output_dir'],
        config['pdb_filter'],
    )

    # MSA
    align(config['output_dir'])

    # Build profile HMM
    build_model(config['output_dir'], config['target']['interpro_id'])

    # Download PDB files for validation
    files_number = download_validation_set(
        config['output_dir'],
        pdb_ids,
        config['validation']['dataset_fraction'],
        config['target']['interpro_id'],
    )

    # Validate
    validate(
        config['output_dir'],
        config['target']['interpro_id'],
        config['target']['interpro_id'],
        config['validation']['e_value_threshold'],
    )

    # Download random PDB files for validation
    download_validation_random_set(config['output_dir'], files_number)

    # "Validate" using random fasta files
    validate(
        config['output_dir'],
        config['target']['interpro_id'],
        "random",
        config['validation']['e_value_threshold'],
    )

    # Calculate metrics
    TP, FP, TN, FN = get_metrics(
        config['output_dir'],
        config['target']['interpro_id'],
        "random",
    )

    # Model didn't pass validation since it is not able to predict.
    if TP == 0:
        print("""
        hmmsearch routine returned 0 correctly predicted structures.
        Something is wrong with the model or config. One of the possible reasons
        - E-value threshold is too string, check your config.
        """)
        sys.exit(1)

    precision= TP/(TP+FP)
    recall = TP/(TP+FN)
    f1 = 2*precision*recall/(precision+recall)

    # Draw a confusion matrix
    matrix = sns.heatmap([
        [TP/(TP+FP+FN+TN), FP/(TP+FP+FN+TN)],
        [FN/(TP+FP+FN+TN), TN/(TP+FP+FN+TN)],
    ], annot=True, fmt='.2%', cmap='Blues')
    matrix.set(xlabel=f"Precision: {round(precision,2)} Recall: {round(recall,2)} F1: {round(f1,2)}")
    fig = matrix.get_figure()
    fig.savefig(f"{config['output_dir']}/validation/confusion_matrix.png")

    # Search
    if config['search']['do']:
        search(
            config['output_dir'],
            config['search']['e_value_threshold'],
            config['target']['interpro_id'],
            config['search']['db_fraction'],
            config['search']['db']
        )
