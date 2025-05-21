!/bin/bash
for img in *.png; do
    filename=${img%.*}
    echo $filename
    convert -trim "$filename.png" "$filename.png"
done


echo **DONE**
