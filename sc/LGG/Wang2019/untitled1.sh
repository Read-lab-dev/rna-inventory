for name in *.gz;
    do 
     echo ${name#*_*_};
done
