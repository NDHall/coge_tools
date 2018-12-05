##GBS Markers and A vs B Bias

###GBS Markers
* Markers taken from Sharma et al. 2018 **Genome wide association mapping of agro-morphological traits among a diverse collection
 of finger millet (*Eleusine coracana* L.) genotypes using SNP markers** 
 
 * These markers were first described in previous publication Kumar et al 2016 **Genotyping-by-Sequencing Analysis for 
 Determining Population Structure of Finger Millet Germplasm of Diverse Origins**
 
 * Markers were downloaded from excel spread sheet on google drive as 2 columns. Marker  name and sequence. 
 
 * awk was used to convert tsv to fasta such that name was fasta header.
 
 * Fasta was mapped against entire genome using `bowtie2 -U -f  --very-fast` on beta. 
 
 * Significant marker by traits were converted into a tsv.
 
 
 