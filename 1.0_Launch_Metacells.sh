# nohup docker exec seacell_aml bash /AML/1.0_Launch_Metacells.sh >/home/sevastopol/data/gserranos/SEACells_AML/AML_seacells_nohup.out 2>&1 &

conda init bash
conda activate seacells
# python -c "import SEACells; SEACells.__version__"
echo "Starting to run the script fot MDS"
# python -u /AML/1.2_Seacells_4_AML_DeNovo.py
python -u /AML/1.0_Launch_Metacells.py