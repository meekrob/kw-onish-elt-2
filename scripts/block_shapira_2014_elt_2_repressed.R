block_shapira_2014_repressed_tx <- c(
 'B0334.2', 'B0350.2', 'B0507.10', 'C02C2.1', 'C05D10.4',
 'C18A11.7', 'C34H4.4', 'C36E6.5', 'C42D8.8', 'C44B11.3',
 'C51E3.7', 'C53C11.3', 'F02G3.1', 'F09F3.6', 'F32B5.6',
 'F38B2.1', 'F44F4.4', 'F49E11.11', 'F53F1.5', 'F54E4.3',
 'F57B1.4', 'K07E1.1', 'K08D8.3', 'K10B3.7', 'R08E3.4',
 'T03F1.11', 'T05H10.3', 'T18D3.4', 'T19A6.2', 'T19C9.8',
 'T21B6.3', 'T22C1.7', 'Y105C5A.12', 'Y37D8A.23', 'Y67A10A.7',
'F17C8.2', 'F32G8.4', 'R06C7.10', 'R12E2.7', 'T19A6.2',
'T22A3.2', 'T22B2.4', 'W05G11.6', 'Y105E8B.1')

paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

repressed_by_elt_2 <-getBM(mart = paramart, 
     filter="wormbase_gseqname",
     value=block_shapira_2014_repressed_tx,
     attributes = c(
        'chromosome_name', 
        'start_position', 
        'end_position', 
        'strand',
        'wormbase_gseq',
        'wbps_gene_id',
        "external_gene_id"))
