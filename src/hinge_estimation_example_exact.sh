echo 'adenylate kinase'

docker exec shilo ./execute_shibuya_method adenylate_kinase/1ake.pdb adenylate_kinase/4ake.pdb A A 4

echo 'short linker loop'
docker exec shilo ./execute_shibuya_method linker/1mdt.pdb linker/1ddt.pdb A A 2

echo 'long linker loop'
docker exec shilo ./execute_shibuya_method linker/1nkh.pdb linker/1fgx.pdb B A 2

echo 'domain swapping'
docker exec shilo ./execute_shibuya_method 2oct.pdb pdb4n6v.pdb A A 3
