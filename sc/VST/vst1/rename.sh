for name in *.gz;
    do 
    	mkdir ${name%%-*}
    	mv $name ${name%%-*}/${name#*_*-};
done
