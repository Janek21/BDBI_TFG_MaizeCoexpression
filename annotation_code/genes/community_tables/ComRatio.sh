

comms=$(grep -v "Community" $1|cut -d$'\t' -f1|uniq) #no sort needed, they should already be ordered

echo "Community"$'\t'"Known"$'\t'"Unknown" > $2
for c in $comms; do
	kw=$(grep $c $1 |grep "Known" --no-ignore-case -c) #known genes
	ukw=$(grep $c $1 |grep "Known" --no-ignore-case -vc) #unknwn genes
	
	echo $c$'\t'$kw$'\t'$ukw >> $2
done
