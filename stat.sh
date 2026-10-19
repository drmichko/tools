echo node 0
ssh rayol  ps -eF  | grep -E '(exe|CMD)' | grep -v root
echo node 1
ssh rayol-node-1    ps -eF  | grep -E '(exe|CMD)'| grep -v root
echo node 2
ssh rayol-node-2    ps -eF  | grep -E '(exe|CMD)'| grep -v root

echo 
ls  -1 rho | tr -cd '[0-9\n]'	| sort  -rn   | tr '\n'  ' '  | fold -s | tail -n 4

echo

