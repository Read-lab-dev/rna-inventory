for name in *.gz;
    do 
    	mkdir ${name:11:2}
    	mv $name ${name:11:2}/${name#*_*_};
done
