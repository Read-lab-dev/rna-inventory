for name in *.gz;
    do 
    	extracted_string=$(echo "$name" | sed 's/^[^_]*_//;s/_.*//')
    	mkdir $extracted_string
    	mv $name $extracted_string/${name#*_*_};
done
