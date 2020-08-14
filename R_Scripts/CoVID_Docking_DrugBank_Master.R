
# Library and Functions ---------------------------------------------------

setwd("~/Documents/covid19-docking-master/")
require(data.table)
library(readxl)
require(ggpubr)
source("~/Documents/covid19-docking-master/scripts/dock_helpers.r")

# Loading data for query --------------------------------------------------

drugbank_query = read_xlsx("2020-05-07  drugbank and pubchem id details.xlsx", sheet = 1)


         ############        DRUG BANK           ###############
# DRUGBANK DOCKING --------------------------------------------------------

df_drugbank = as.data.frame(drugbank_query[,c("drugbank_id","Name in drugbank","smiles")])
colnames(df_drugbank) = c('drugbank_id', 'name', 'smiles')
#Check if any N.A values are present
any(is.na(df_drugbank))


drugbank_dock = function(i){
  message(i)
  name = df_drugbank$drugbank_id[i]
  smiles = df_drugbank$smiles[i]
  
  path = generate_ligand_pdbqt(smiles, name, out_dir = '~/Documents/DrugBankQueryLigands/pdbqt')
  
  if(is.na(path)){
    unlink(sprintf('~/Documents/DrugBankQueryLigands/pdbqt/%s.pdbqt', name))
    unlink(sprintf('~/Documents/DrugBankQueryLigands/docked/docked_%s.pdbqt', name))
    unlink(sprintf('~/Documents/DrugBankQueryLigands/logs/%s.pdbqt.log', name))
    return(NULL)
  }
  
  fout = sprintf('~/Documents/DrugBankQueryLigands/docked/docked_%s', basename(path))
  if(file.exists(fout)){
    message('-- skipping ', name)
    return(NULL)
  }
  
  run_docking( path, exhaustiveness = 10, cores = 1, active_site = T, 
               dock_dir = '~/Documents/DrugBankQueryLigands/docked', log_dir = '~/Documents/DrugBankQueryLigands/logs' )
}



# Running Docking for every 7 ligands -------------------------------------

lapply(1:7, drugbank_dock)
lapply(7:14, drugbank_dock)
lapply(14:21, drugbank_dock)
lapply(21:28, drugbank_dock)


# Running zero coordinate for pdbs ----------------------------------------

source("~/Documents/covid19-docking-master/scripts/find_zero_coordinate_pdb.r")


# Making md table ---------------------------------------------------------

setwd("~/Documents/DrugBankQueryLigands/")
x = list.files('logs/', full.names = T)


df = rbindlist(lapply(x, function(file){
  z = readLines(file)
  sp = strsplit(grep('^1 ', trimws(z), value = T), '\\s+')
  if(length(sp) == 0) return(NULL)
  v = sp[[1]][2]
  v = as.numeric(v)
  data.table(file = basename(file), value=v)
}))

df$drugbank_id = gsub('\\..*', '', df$file)

x = df_drugbank[,c('drugbank_id', 'name')]
m = merge(df, x, by='drugbank_id')
m = m[order(value)]

q = subset(m, name == 'Quercetin')
p = gghistogram(subset(m, value <0), 'value', bins = 50, fill='#333333', color='white') + grids(color='#eeeeee') + 
  xlab('Affinity kcal/mol') + ylab('Count') + geom_vline(xintercept = q$value, linetype=2, color='red') + 
  ggtitle('Distribution of affinity values for DrugBank subset compounds')
#print(p)


tab = m[,c('name', 'drugbank_id', 'value')]

# Write table
tmp = tab
names(tmp) = c('Compound', 'DrugBank ID', 'Affinity')
write.csv(tmp, 'drugbank_docking_affinity_results.csv', row.names=F)

tab$drugbank_id = sprintf('[%s](https://www.drugbank.ca/drugs/%s)', tab$drugbank_id, tab$drugbank_id)
ind = nchar(tab$name) > 80
tab$name[ind] = paste0(substring(tab$name[ind], 1, 77), '...')
names(tab) = c('Compound', 'DrugBank ID', 'Affinity')


ggsave('images/hist_drugbank.png', p)

writeLines(knitr::kable(head(tab, 20)), 'images/top20.md')


# Generating pymol images -------------------------------------------------


run_pymol <- function(ligand_path, image_path){
  
  script = sprintf('
reinitialize 

##### Load protein
load ~/Documents/DrugBankQueryLigands/protein/protein_6yb7.pdbqt, protein
select hetatm
remove sele
show surface, protein

# Hide cartoon and set surface
hide cartoon
set surface_color, white
bg_color white
#set solvent_radius, 1

load %s, drug
show spheres, drug
set sphere_scale, 0.7, drug

set_view (\\
     0.653096735,    0.308816165,   -0.691443324,\\
    -0.395566553,    0.917719483,    0.036244851,\\
     0.645744741,    0.249840796,    0.721521676,\\
    -0.000025807,    0.000011973, -173.567214966,\\
     2.680393219,   -2.178984165,   22.382369995,\\
   115.596817017,  231.541671753,  -20.000000000 )
   
# Color the drug
color rutherfordium, elem c
color platinum, elem h
color meitnerium, elem o
color rubidium, elem n
color orange, elem s

##### Ray trace options
set ray_trace_mode, 1
set ray_trace_fog, 0
# set ray_shadows, 0
unset depth_cue
set antialias,2

ray 1600,1200
png %s, 1600, 1200

', ligand_path, image_path)
  
  writeLines(script, 'script.pml')
  system('/Users/pradeep/anaconda3/bin/pymol -cq script.pml')
  unlink('script.pml')
  
}

# pymol drugbank
db = head(fread('~/Documents/DrugBankQueryLigands/drugbank_docking_affinity_results.csv'), 20)
picked_db = sprintf('~/Documents/DrugBankQueryLigands/docked/docked_%s.pdbqt', db$`DrugBank ID`)
picked_db_out = sprintf('~/Documents/DrugBankQueryLigands/images/%s.png', basename(picked_db))

for(i in 1:length(picked_db)){
  run_pymol(picked_db[i], picked_db_out[i])
}

