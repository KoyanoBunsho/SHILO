def convert_pdb_id(pdb_str):
    pdb, chain = pdb_str.split("_")
    return f"{pdb.lower()}_{chain}"
