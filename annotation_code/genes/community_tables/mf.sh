

allC=$(grep -v 'Community' turquoise_commF.txt|cut -f1|sort|uniq)

for c in $allC; do
	echo $c
	grep $c turquoise_commF.txt |cut -f4 -d$'\t'|cut -d'.' -f1|sort|uniq -c|sort -n|tail -n9|grep 'e'|grep -v '_'
done

