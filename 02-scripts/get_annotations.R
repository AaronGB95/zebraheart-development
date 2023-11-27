

# Paquetes necesarios
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Dr.eg.db)


# Paquetes org.XX.eg.db

## Rutas KEGG
keytypes(org.Dr.eg.db) # identificadores disponibles en el paquete

# Por defecto el identificador que usa es el Entrezid
# toTable da la información en forma de dataframe

id2name <- toTable(org.Dr.egGENENAME) 

# Obtener las rutas KEGG de cada gen con el identificador Gene symbol

## 1. Lista de genes 
keys <- keys(org.Dr.eg.db, keytype="SYMBOL")

## 2. Asociar gen con pathway
### Argumentos
### keys -> genes de interés
### columns -> identifadores de interés
### keytype -> tipo de identificador de tus genes de interés
r_kegg <- AnnotationDbi::select(org.Dr.eg.db, keys=keys,
                                        columns=c("SYMBOL", "PATH"), 
                                        keytype="SYMBOL")

## 3. Eliminar genes sin anotación
r_kegg <- na.exclude(r_kegg) 

## 4. Comprobaciones
table(duplicated(r_kegg)) 
head(r_kegg)

# Obtener rutas KEGG con los identificadores de Ensembl
keys <- keys(org.Dr.eg.db, keytype="ENSEMBL")
rat_ensembl_kegg <- AnnotationDbi::select(org.Dr.eg.db, keys= keys,
                                        columns=c("ENSEMBL", "PATH"), keytype="ENSEMBL")
rat_ensembl_kegg <- na.exclude(rat_ensembl_kegg) 
table(duplicated(rat_ensembl_kegg))
head(rat_ensembl_kegg)


## Términos GO

keys <- keys(org.Dr.eg.db, keytype="SYMBOL")
rat_go <- AnnotationDbi::select(org.Dr.eg.db, keys= keys,
                               columns=c("SYMBOL", "GO"), keytype="SYMBOL")

rat_go <- na.exclude(rat_go)

# No solo se incluye el identificador del término GO, también incluye el 
# código de evidencia y la ontología.

# Filtrar por ontologías

## Biological process
r_bp <- rat_go[rat_go$ONTOLOGY == "BP",]
r_bp <- r_bp[,1:2]
table(duplicated(r_bp))
r_bp <- r_bp[!duplicated(r_bp), ]

## Biological process
r_mf <- rat_go[rat_go$ONTOLOGY == "MF",]
r_mf <- r_mf[,1:2]
table(duplicated(r_mf))
r_mf <- r_mf[!duplicated(r_mf), ]

## Cellular component
r_cc <- rat_go[rat_go$ONTOLOGY == "CC",]
r_cc <- r_cc[,1:2]
table(duplicated(r_cc))
r_cc <- r_cc[!duplicated(r_cc), ]


# Paquete biomaRt
# A veces no puede conectarse al servidor y da error. Si pasa esto
# hay que cambiar de host:
# "uswest.ensembl.org", "asia.ensembl.org", "useast.ensembl.org", "www.ensembl.org",

mart <- useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl")

# mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl",
#                host="uswest.ensembl.org")

# Ver atributos y filtros disponibles
View(as.data.frame(listAttributes(mart)))
View(as.data.frame(listFilters(mart)))

# Obtener rutas de reactome de todos los genes
r_reac <- getBM(attributes=c("go_id", "name_1006", "namespace_1003"), 
               values = TRUE , mart = mart)

# Eliminar genes sin anotación
r_reac <- r_reac[!(r_reac$reactome == ""),  ]
table(duplicated(r_reac))


# Escribir resultados

# Lista de resultados
list_results <- list(r_bp, r_mf, r_cc, r_reac, r_kegg)
names(list_results) <- c("rat_symbol_GO.BP", "rat_symbol_GO.MF", "rat_symbol_GO.CC", 
                        "rat_symbol_Reactome", "rat_symbol_KEGG")

# Crear directorio de anotaciones
dir.create("./annotations")

# Escribir resutlados
lapply(1:length(list_results), function(i){
  file <- paste0("./annotations/", names(list_results)[i], ".txt")
  write.table(list_results[[i]], file, sep = "\t", quote = F,  row.names = F, col.names = T)
})

