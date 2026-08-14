1. bash extractUnikmers.sh prefix 21
2. nohup python -u read_unikmer_map.py > read_unikmer_map.log 2>&1 &
3. nohup python -u bin_reads.py > bin_reads.log 2>&1 &
4. nohup python -u universe.py > universe.log 2>&1 &
5. nohup python -u anchor.py > anchor.log 2>&1 &
6. nohup python -u direct_bridge.py > direct_bridge.log 2>&1 &
7. nohup python -u indexer.py > indexer.log 2>&1 &
8. nohup python -u barcode_assembler.py > barcode_assembler.log 2>&1 &
9. nohup python -u final_read_selection.py > final_read_selection.log 2>&1 &
10. bash assembly.sh

cd ..

11. quast.py /24-2/home/vselv001/MTPLite/output/mtpv1.1/mtpv1.1.asm.fasta \
  -r /24-2/home/vselv001/MTPLite/reference/chr1.fasta \
  -o /24-2/home/vselv001/MTPLite/output/quast_mtpv1.1

So far @mtp_lite_v1.1 is headed towards the correct direction. The less contigs the better.