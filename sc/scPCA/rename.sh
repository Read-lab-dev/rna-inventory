for name in *.gz;
    do 
    	mkdir ${name%%_*}
    	mv $name ${name%%_*}/${name#*_*_};
done
