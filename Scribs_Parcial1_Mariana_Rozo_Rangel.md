Parcial.

##Primer punto
scp -i bio.pt.pem -P 53841 bio.pt@loginpub-hpc.urosario.edu.co:/home/bio.pt/data/marianarozor/clases/ensamblaje/Pacbio/quast_hifi_ensamblaje/ /mnt/c/Users/Mariana\ Rozo/Downloads/
(base) [bio.pt@caldas01 datos]$ cp -r SRR11140744/ /home/bio.pt/data/Parcial/parcial1_marianarozorangel/
(base) [bio.pt@caldas01 datos]$ cp *.fasta /home/bio.pt/data/Parcial/parcial1_marianarozorangel/

#PASAR DE FORMATO .sra a fastq
#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1 # Numero de nodos
#SBATCH -n 4 # Numero de nucleos
#SBATCH -t 0-00:20 # limite de tiempo (D-HH:MM)
#SBATCH -o salida1.out # Salida STDOUT
#SBATCH -e error1.err # Salida STDERR
#SBATCH --mail-user=mariana.rozor@urosario.edu.co #Direccion e-mail a donde notificar el estado del trabajo
#SBATCH --mail-type=ALL #Especifica que eventos notificar al correo (ALL, BEGIN, END, REQUEUE, FAIL)
#SBATCH --reservation=reserva_dev
module load sratoolkit/3.2.0
fasterq-dump SRR11140744/SRR11140744.sra --split-files


### FASTQC
#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1 # Numero de nodos
#SBATCH -n 4 # Numero de nucleos
#SBATCH -t 0-00:20 # limite de tiempo (D-HH:MM)
#SBATCH -o salida1.out # Salida STDOUT
#SBATCH -e error1.err # Salida STDERR
#SBATCH --mail-user=mariana.rozor@urosario.edu.co #Direccion e-mail a donde notificar el estado del trabajo
#SBATCH --mail-type=ALL #Especifica que eventos notificar al correo (ALL, BEGIN, END, REQUEUE, FAIL)
#SBATCH --reservation=reserva_dev
module load fastqc
fastqc *.fastq


scp -i bio.pt.pem -P 53841 bio.pt@loginpub-hpc.urosario.edu.co:/home/bio.pt/data/Parcial/parcial1_marianarozorangel/SRR11140744_1_fastqc.html /mnt/c/Users/Mariana\ Rozo/Downloads/
SRR11140744_1_fastqc.html                                                             100%  241KB 777.7KB/s   00:00




#### Kmers antes de filtrado

#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1 # Numero de nodos
#SBATCH -n 4 # Numero de nucleos
#SBATCH -t 0-00:20 # limite de tiempo (D-HH:MM)
#SBATCH -o salida1.out # Salida STDOUT
#SBATCH -e error1.err # Salida STDERR
#SBATCH --mail-user=mariana.rozor@urosario.edu.co #Direccion e-mail a donde notificar el estado del trabajo
#SBATCH --mail-type=ALL #Especifica que eventos notificar al correo (ALL, BEGIN, END, REQUEUE, FAIL)
#SBATCH --reservation=reserva_dev

module load jellyfish.

jellyfish count -C -m 21 -s 100M -t 4 \
SRR11140744_1.fastq SRR11140744_2.fastq \
-o before.jf

jellyfish histo before.jf > beforekmer1.histo


## DESPUES KMER.

#!/bin/bash
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-00:20
#SBATCH -o salidakmer_after.out
#SBATCH -e errorkmer_after.err
#SBATCH --mail-user=mariana.rozor@urosario.edu.co
#SBATCH --mail-type=ALL
#SBATCH --reservation=reserva_dev

module load java11

java -jar /opt/ohpc/pub/apps/trimmomatic/0.38/bin/trimmomatic-0.38.jar PE \
-phred33 \
SRR11140744_1.fastq SRR11140744_2.fastq \
R1_paired.fastq R1_unpaired.fastq \
R2_paired.fastq R2_unpaired.fastq \
ILLUMINACLIP:/opt/ohpc/pub/apps/trimmomatic/0.38/bin/adapters/TruSeq3-PE.fa:2:30:10
LEADING:3 TRAILING:3
SLIDINGWINDOW:4:20
MINLEN:50

#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-00:20
#SBATCH -o salidakmer_after.out
#SBATCH -e errorkmer_after.err
#SBATCH --mail-user=mariana.rozor@urosario.edu.co
#SBATCH --mail-type=ALL
#SBATCH --reservation=reserva_dev

module load spades

spades.py --careful -1 R1_paired.fastq -2 R2_paired.fastq -o SARSCOV_illumina



/datacnmat01/ciencias/appsbio/conda/pkgs/quast-5.2.0-py39pl5321h4e691d4_3/opt/quast-5.2.0/quast.py scaffolds.fasta -o quast_results



Punto 2

#CONCATENADO:

#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-00:20
#SBATCH -o salida.out
#SBATCH -e error.err
#SBATCH --mail-user=mariana.rozor@urosario.edu.co
#SBATCH --mail-type=ALL
#SBATCH --reservation=reserva_dev

cat *.fasta >> allspin_unido.fasta


#Creacion de base y comparacion:

#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-00:20
#SBATCH -o salida.out
#SBATCH -e error.err
#SBATCH --mail-user=mariana.rozor@urosario.edu.co
#SBATCH --mail-type=ALL

module load muscle/3.8.31
muscle -in allspin_unido.fasta -out aln_sequences_muscle.fasta


#PUNTO 3
wget "https://raw.githubusercontent.com/paula-torres/Bioinfo_UR/refs/heads/main/files/df_master.csv"

library(ggplot2)
#MEJORAR GRAFICOS
datos<-read.csv('df_master.csv.1',head=T)
mynamestheme <- theme(plot.title = element_text(size = 16, hjust = 0.5))

##Una gráfica de puntos con la Fecha en el eje X y las Mutaciones en el eje Y, estos puntos deben estar clasificados por colores ##según la cepa.

#(base) [bio.pt@caldas01 tablapunto3]$ head df_master.csv
#,Fecha,t_num,Clado,Pico_Frec,Frecuencia,Casos_Totales,Mutaciones
#1,2019-12-01,18231,Alpha/Beta/Gamma,18659,3.65E-05,16.85479224,NA
#scp -i bio.pt.pem -P 53841 bio.pt@loginpub-hpc.urosario.edu.co:/home/bio.pt/data/Parcial/parcial1_marianarozorangel/tablapunto3/df_master.csv /mnt/c/Users/Mariana\ Rozo/Downloads/


pdf("Grafico1.pdf")

plot_20 <- ggplot(datos, aes(x =Fecha , y = Mutaciones, color = Clado )) +
  geom_point() +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  mynamestheme +
  ggtitle("Variación de las mutaciones a traves del tiempo") +
  ylab("Fecha") +
  xlab("Numero de Mutaciones")


# Una gráfica de líneas que relacione la Fecha (eje X) con el número de casos (eje Y).
pdf("Grafico2.pdf")
plot_21 <- ggplot(datos, aes(x = Fecha, y = Casos_Totales , color = Clado)) +
    geom_line() +
    theme(panel.background = element_rect(fill = "white", colour = "black")) +
    mynamestheme +
    ggtitle("Numero de casos por fecha") +
    ylab("Numero de casos") +
    xlab("Fecha")
print(plot_20)
dev.off()
