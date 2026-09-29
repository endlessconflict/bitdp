# References

Work this library builds on or is compared against. The code was written from the mathematics in these sources, not from any existing implementation.

## Papers

- G. Myers. A fast bit-vector algorithm for approximate string matching based on dynamic programming. *Journal of the ACM* 46(3), 1999. doi:10.1145/316542.316550. The unit-cost kernel that `reference.myers` implements and that the derived kernel is compared against.
- J. Loving, Y. Hernandez, G. Benson. BitPAl: a bit-parallel, general integer-scoring sequence alignment algorithm. *Bioinformatics* 30(22):3166–3173, 2014. doi:10.1093/bioinformatics/btu507. Prior bit-parallel work for match/mismatch/gap weights: a construction designed by hand for that family (one word per difference value, zones split at I - G, additions through runs of the minimum difference), generated per weight set.
- T. Mytkowicz, M. Musuvathi, W. Schulte. Data-parallel finite-state machines. *ASPLOS* 2014. doi:10.1145/2541940.2541988. Parallel evaluation of finite-state machines by composing their transition functions, the idea behind reading a DP column as a machine.

- K. E. Batcher. Sorting networks and their applications. *AFIPS Spring Joint Computer Conference* 32:307–314, 1968. doi:10.1145/1468075.1468121. The odd-even merging network used for unary addition.

## Books

- G. Gopalakrishnan. *Computation Engineering: Applied Automata Theory and Logic*. Springer, 2006. doi:10.1007/0-387-32520-4. Finite-state machines, transducers and boolean function representation.
- D. E. Knuth. *The Art of Computer Programming, Volume 4, Fascicle 6: Satisfiability*. Addison-Wesley, 2015. Cardinality constraints encoded with sorting networks, the source of the idea that merging sorted bit strings counts in unary.
