<p align="center">
<img src="../logo.svg" alt="DNA Spider-Web" title="DNASpiderWeb" width="40%"/>
</p>

---

You can generate the results (25MB) by yourself (step 2 in pipeline folder), it only takes several hours.


Ignored file list:
```html
├── VLC1[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 1
├── VLC2[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 2
├── VLC3[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 3
├── VLC4[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 4
├── VLC5[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 5
├── VLC6[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 6
├── VLC7[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 7
├── VLC8[t].npy     // Transcoding results of variable-length graph coder with biochemical filter 8
```

The local biochemical constraint groups are as follow:

|group index | maximum homopolymer runs | GC content range| undesired DNA motifs |
| ---- | ---- | ---- | ---- |
| 1 |  2  | 50%~50% |  AGCT, GACGC, CAGCAG, GATATC, GGTACC, CTGCAG, GAGCTC, GTCGAC, AGTACT, ACTAGT, GCATGC, AGGCCT, TCTAGA |
| 2 |  1  |   N/A   | N/A |
| 3 | N/A | 10%~30% | N/A |
| 4 |  2  | 40%~60% | AGA, GAG, CTC, TCT |
| 5 |  2  | 40%~60% | N/A |
| 6 | N/A | 50%~70% | N/A |
| 7 |  4  | 40%~60% | N/A |
| 8 |  3  |   N/A   | N/A |

The results of fixed-length code are as follows:

<p align="center">
<img src="fixed_results.png" alt="results of fixed-length code" title="generated fixed" width="100%"/>
</p>


And those of variable-length code are as follows:

<p align="center">
<img src="variable_results.png" alt="results of variable-length code" title="generated variable" width="100%"/>
</p>