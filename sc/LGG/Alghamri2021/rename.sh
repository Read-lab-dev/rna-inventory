for name in *.gz;
    do 
    	extracted_string=${name%%_*}
    	mkdir $extracted_string
    	mv $name $extracted_string/${name#*_*_*_*_};
done
