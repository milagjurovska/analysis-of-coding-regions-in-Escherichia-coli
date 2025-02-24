install.packages("rentrez")
library(rentrez)

genbank_id <- "NC_000913"
seq <- entrez_fetch(db = "nucleotide", id = genbank_id, rettype = "fasta", retmode = "text")

cat(seq)

fasta_sequence <- strsplit(seq, "\n")[[1]]
dna_sequence <- paste(fasta_sequence[2:length(fasta_sequence)], collapse = "")
short_sequence <- substr(dna_sequence, 1, 1000)

numeric_sequence <- chartr("ATCG", "1234", short_sequence)

numeric_sequence <- as.numeric(strsplit(numeric_sequence, "")[[1]])

fft_result <- fft(numeric_sequence)
fft_result

n <- length(numeric_sequence)
freqs <- (0:(n - 1)) / n

target_frequency <- 1 / 3
target_index <- which.min(abs(freqs - target_frequency))

plot(freqs, Mod(fft_result), type = "h",
     main = "Фуриева анализа на ДНК секвенца",
     xlab = "Фреквенција", ylab = "Големина",
     ylim = c(0,100))

abline(v = target_frequency, col = "red", lwd = 2, lty = 2)
text(target_frequency, max(Mod(fft_result)), labels = "f = 1/3", pos = 4, col = "red")
