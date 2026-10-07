#!/bin/bash
echo "====================================================================="
echo " SCRIPT 1: AUTO-POCKET DETECT + UNIFIED PDF REPORT"
echo "====================================================================="
mkdir -p 01_raw_data 02_protein_prep 03_ligand_prep 04_docking 05_results 06_visualization client_delivery

echo "==> Batch Converting All Ligands..."
shopt -s nullglob
for ligand_file in 01_raw_data/*.{sdf,mol2,pdb}; do
    if [ -f "$ligand_file" ]; then
        fname_base=$(basename "$ligand_file")
        fname="${fname_base%.*}"
        obabel "$ligand_file" -O "03_ligand_prep/${fname}.pdbqt" -h --gen3D 2>/dev/null
    fi
done

for receptor_path in 01_raw_data/*.pdb; do
    [ ! -f "$receptor_path" ] && continue
    RECEPTOR_NAME=$(basename "$receptor_path")
    prot_base=$(basename "$RECEPTOR_NAME" .pdb)
    echo "========================================================"
    echo " PROCESSING RECEPTOR: ${RECEPTOR_NAME}"
    echo "========================================================"
    rm -rf 02_protein_prep/* 04_docking/* 05_results/* 06_visualization/*
    
    echo "==> [1] Cleaning PDB (Extracting pure ATOMs only)..."
    grep "^ATOM" "$receptor_path" > "02_protein_prep/${RECEPTOR_NAME}"
    
    echo "==> [2] Running fpocket for automated cavity detection..."
    cd 02_protein_prep
    fpocket -f "${RECEPTOR_NAME}" > /dev/null 2>&1
    out_pdb="${prot_base}_out/${prot_base}_out.pdb"
    coords=$(python3 -c '
import sys
try:
    xs, ys, zs = [], [], []
    with open(sys.argv[1]) as f:
        for line in f:
            if "STP" in line and line[22:26].strip() == "1":
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    if xs: print(f"{sum(xs)/len(xs):.3f} {sum(ys)/len(ys):.3f} {sum(zs)/len(zs):.3f}")
    else: print("0.000 0.000 0.000")
except: print("0.000 0.000 0.000")' "$out_pdb")
    center_x=$(echo "$coords" | awk '{print $1}')
    center_y=$(echo "$coords" | awk '{print $2}')
    center_z=$(echo "$coords" | awk '{print $3}')
    cd ..
    
    if [ "$center_x" = "0.000" ] && [ "$center_y" = "0.000" ]; then
        echo "    -> ERROR: fpocket failed to detect a binding site. Skipping ${RECEPTOR_NAME}..."
        continue
    fi
    echo "    -> Found Pocket: X=${center_x}, Y=${center_y}, Z=${center_z}"
    
    echo "==> [3] Preparing Receptor PDBQT..."
    cd 02_protein_prep
    obabel "${RECEPTOR_NAME}" -O "receptor.pdbqt" -xr -h 2>/dev/null
    obabel "${RECEPTOR_NAME}" -O "receptor_clean.pdb" -xr 2>/dev/null
    cd ..
    
    cat << INN > 04_docking/conf.txt
receptor = 02_protein_prep/receptor.pdbqt
center_x = ${center_x}
center_y = ${center_y}
center_z = ${center_z}
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 8
INN

    echo "==> [4] Executing Batch Docking..."
    for lig_pdbqt in 03_ligand_prep/*.pdbqt; do
        if [ -f "$lig_pdbqt" ]; then
            lig_name=$(basename "$lig_pdbqt" .pdbqt)
            vina --config 04_docking/conf.txt --ligand "$lig_pdbqt" --out "05_results/${lig_name}_posed.pdbqt" > "04_docking/${lig_name}_vina" 2>&1
        fi
    done
    
    echo "==> [5] Running PLIP Profiling & Compiling Pristine 2-Page PDF..."
    for posed_file in 05_results/*_posed.pdbqt; do
        if [ -f "$posed_file" ]; then
            pname=$(basename "$posed_file" _posed.pdbqt)
            # CRITICAL FIX: Extract only the top pose to prevent PLIP stacking crash
            obabel "$posed_file" -O "05_results/${pname}_posed.pdb" -f 1 -l 1 2>/dev/null
            grep -v "^END" 02_protein_prep/receptor_clean.pdb > "06_visualization/${pname}_complex.pdb"
            grep -v "^END" "05_results/${pname}_posed.pdb" >> "06_visualization/${pname}_complex.pdb"
            echo "END" >> "06_visualization/${pname}_complex.pdb"
            mkdir -p "06_visualization/${pname}_tmp"
            
            plip -f "06_visualization/${pname}_complex.pdb" -y -t -o "06_visualization/${pname}_tmp" >/dev/null 2>&1
            find "06_visualization/${pname}_tmp" -name "*.txt" -exec cp {} "06_visualization/${pname}_report.txt" \; 2>/dev/null || true
            find "06_visualization/${pname}_tmp" -name "*.pse" -exec cp {} "06_visualization/${pname}.pse" \; 2>/dev/null || true
            rm -rf "06_visualization/${pname}_tmp"
            
            # CRITICAL FIX: Secure sys.argv parameter passing for python script
            python3 -c "
import sys, os
try: from fpdf import FPDF
except ImportError: sys.exit(0)
pname = sys.argv[1]
class PDF(FPDF):
    def footer(self):
        self.set_y(-10)
        self.set_font('Helvetica', 'I', 8)
        self.cell(0, 5, f'Page {self.page_no()}', 0, 0, 'C')
try:
    pdf = PDF(orientation='L', unit='mm', format='A4')
    pdf.set_margins(left=8, top=10, right=8)
    pdf.add_page()
    pdf.set_font('Helvetica', 'B', 10)
    pdf.cell(0, 6, f'Molecular Docking & Interaction Report: {pname}', 0, 1, 'C')
    pdf.ln(2)
    pdf.set_font('Helvetica', 'B', 9)
    pdf.cell(0, 5, '1. AutoDock Vina Scoring & Log', 0, 1, 'L')
    pdf.set_font('Courier', '', 5.5)
    vina_path = f'04_docking/{pname}_vina'
    if os.path.exists(vina_path):
        with open(vina_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                cl = line.rstrip().encode('latin-1', 'replace').decode('latin-1')
                if len(cl) > 160: cl = cl[:160]
                pdf.cell(0, 3.3, cl, ln=1)
    pdf.add_page()
    pdf.set_font('Helvetica', 'B', 10)
    pdf.cell(0, 6, f'Molecular Docking & Interaction Report: {pname}', 0, 1, 'C')
    pdf.ln(2)
    pdf.set_font('Helvetica', 'B', 9)
    pdf.cell(0, 5, '2. PLIP Interaction Profile', 0, 1, 'L')
    pdf.set_font('Courier', '', 5.5)
    plip_path = f'06_visualization/{pname}_report.txt'
    if os.path.exists(plip_path):
        with open(plip_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                cl = line.rstrip().encode('latin-1', 'replace').decode('latin-1')
                if len(cl) > 160: cl = cl[:160]
                pdf.cell(0, 3.3, cl, ln=1)
    pdf.output(f'06_visualization/{pname}_report.pdf')
except: pass" "$pname"
        fi
    done
    
    echo "==> [6] Organizing final delivery package..."
    mkdir -p "client_delivery/${prot_base}_delivery"
    cp 05_results/*_posed.pdbqt "client_delivery/${prot_base}_delivery/" 2>/dev/null || true
    cp 04_docking/*_vina "client_delivery/${prot_base}_delivery/" 2>/dev/null || true
    cp 06_visualization/*.pse "client_delivery/${prot_base}_delivery/" 2>/dev/null || true
    cp 06_visualization/*_report.txt "client_delivery/${prot_base}_delivery/" 2>/dev/null || true
    cp 06_visualization/*_report.pdf "client_delivery/${prot_base}_delivery/" 2>/dev/null || true
    echo "Pipeline Complete for ${RECEPTOR_NAME}."
done
