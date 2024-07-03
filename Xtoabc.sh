tr -d '* '  < $1 > /tmp/file.txt
sed -i 's/x10/j/g' /tmp/file.txt
sed -i 's/x9/i/g' /tmp/file.txt
sed -i 's/x8/h/g' /tmp/file.txt
sed -i 's/x7/g/g' /tmp/file.txt
sed -i 's/x6/f/g' /tmp/file.txt
sed -i 's/x5/e/g' /tmp/file.txt
sed -i 's/x4/d/g' /tmp/file.txt
sed -i 's/x3/c/g' /tmp/file.txt
sed -i 's/x2/b/g' /tmp/file.txt
sed -i 's/x1/a/g' /tmp/file.txt
sed -i 's/^/anf=/' /tmp/file.txt
cat /tmp/file.txt
