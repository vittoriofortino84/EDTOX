#!/bin/bash
# retireving the tags for KEGG pathways
rm -rf *.tag
rm -rf pathway.html 
wget https://www.genome.jp/kegg/pathway.html
file=pathway.html
cat $file|sed '/class/d'|grep -n \<h4|awk -F . {'print $2'}|sed -e 's/>//g' -e 's/<//g' -e 's/\///g' -e 's/h4//g' >category_name.txtt
cat $file|grep -n \<h4|sed '/class/d'|awk -F : {'print $1'} >beginnings.txtt
n_categories=`cat beginnings.txtt |wc -l`
nn_c=`echo $n_categories -1|bc` 
cat beginnings.txtt |tail -n $nn_c >endings.txtt
cat pathway.html|wc -l >>endings.txtt

for i in `seq 1 $n_categories `
do
cat_name=`cat category_name.txtt|head -n $i|tail -n 1`
fn=`echo $cat_name|sed -e "s/ /_/g"`
echo $fn
start_of_line=`cat beginnings.txtt|head -n $i|tail -n 1`
end_of_line=`cat endings.txtt|head -n $i|tail -n 1`
delta=`echo $end_of_line - $start_of_line|bc`
tt=`cat $file |head -n $end_of_line|tail -n $delta|grep dt|grep show_pathway|sed '/h=/d'|awk -F dt {'print $2'}|sed -e 's/>//g' -e 's/<//g' -e "s/\///g"` 
echo $tt >$fn.tag
#### For subcategorization
cat $file |head -n $end_of_line|tail -n $delta >new.html 
cat new.html |grep -n '<b>'|awk -F : {'print $1'} >new_beginning.txtt
cat new.html |grep -n \<b\>|awk  {'print $2$3$4$5'}|sed -e 's/>//g' -e 's/<//g' -e 's/\/b//g' -e 's/:/_/g' >new_category_name.txtt
n_new_categories=`cat new_beginning.txtt |wc -l`
new_nn_c=`echo $n_new_categories -1|bc`
cat new_beginning.txtt |tail -n $new_nn_c >new_endings.txtt
cat new.html|wc -l >>new_endings.txtt
mkdir -p $fn
rm -rf ./$fn/*.tag
      for j in `seq 1 $n_new_categories `
      do
      new_cat_name=`cat new_category_name.txtt|head -n $j|tail -n 1`
      new_fn=`echo $new_cat_name|sed -e "s/ /_/g"`
      echo $new_fn
      new_start_of_line=`cat new_beginning.txtt|head -n $j|tail -n 1`
      new_end_of_line=`cat new_endings.txtt|head -n $j|tail -n 1`
      new_delta=`echo $new_end_of_line - $new_start_of_line|bc`
      
      new_tt=`cat new.html |head -n $new_end_of_line|tail -n $new_delta|grep dt|grep show_pathway|sed '/h=/d'|awk -F dt {'print $2'}|sed -e 's/>//g' -e 's/<//g' -e "s/\///g"`
      echo $new_tt >./$fn/$new_fn.tag 
      done

done
rm -rf *.txtt 
rm -rf *.html

