for name in *.gz;
    do 
    	mkdir ${name:11:3}
    	mv $name ${name:11:3}/${name#*_*_};
done
