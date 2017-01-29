for i in {1..40}; do
	mct=$(printf %03d "$i")	
	./getchaindynamics.sh $mct
done
