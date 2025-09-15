import json

def convertjson_from_cellvit_to_scMTOP(path_to_cellvit_json,
                                       path_to_scMTOP_json):
    with open(path_to_cellvit_json, 'r') as f:
        cellvit_json = json.load(f)

    mag = cellvit_json['wsi_metadata']['magnification']
    nuc = {}
    for i, item in enumerate(cellvit_json['cells']):
        nuc[i] = item
    
    scMTOP_json = {'mag': mag, 'nuc': nuc}

    with open(path_to_scMTOP_json, 'w') as f:
        json.dump(scMTOP_json, f)
    
