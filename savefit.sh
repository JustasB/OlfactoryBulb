branch_name=fit-`date +%Y%m%d%H%M%S`
git branch $branch_name
git checkout $branch_name

python fit2.py &> fit.log

git add .
git commit -m 'fit results'
git push

