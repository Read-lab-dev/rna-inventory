for name in *.gz;
    do 
    	mkdir ${name:0:10}
    	mv $name ${name:0:10}/${name#*_*_};
done
