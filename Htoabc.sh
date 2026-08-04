tr -d '* '  < $1 > /tmp/file.txt
sed -i 's/9/j/g' /tmp/file.txt
sed -i 's/8/i/g' /tmp/file.txt
sed -i 's/7/h/g' /tmp/file.txt
sed -i 's/6/g/g' /tmp/file.txt
sed -i 's/5/f/g' /tmp/file.txt
sed -i 's/4/e/g' /tmp/file.txt
sed -i 's/3/d/g' /tmp/file.txt
sed -i 's/2/c/g' /tmp/file.txt
sed -i 's/1/b/g' /tmp/file.txt
sed -i 's/0/a/g' /tmp/file.txt
sed -i 's/^/anf=/' /tmp/file.txt
cat /tmp/file.txt
