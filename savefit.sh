branch_name=fit-`date +%Y%m%d%H%M%S`
git branch $branch_name
git checkout $branch_name

jupyter nbconvert --execute --to notebook fitting.ipynb

git add .
git commit -m 'fit results'

#git remote set-url origin https://justasb:THaswu98@github.com/JustasB/OlfactoryBulb.git
#git push origin $branch_name
#shutdown -h
