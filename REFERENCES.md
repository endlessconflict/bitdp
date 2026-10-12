# References

Work this library builds on or is compared against. The code was written from the mathematics in these sources, not from any existing implementation.

## Papers

- G. Myers. A fast bit-vector algorithm for approximate string matching based on dynamic programming. *Journal of the ACM* 46(3), 1999. doi:10.1145/316542.316550. The unit-cost kernel that `reference.myers` implements and that the derived kernel is compared against.
- J. Loving, Y. Hernandez, G. Benson. BitPAl: a bit-parallel, general integer-scoring sequence alignment algorithm. *Bioinformatics* 30(22):3166–3173, 2014. doi:10.1093/bioinformatics/btu507. Prior bit-parallel work for match/mismatch/gap weights: a construction designed by hand for that family (one word per difference value, zones split at I - G, additions through runs of the minimum difference), generated per weight set.
- A. Bergeron, S. Hamel. Vector algorithms for approximate string matching. *International Journal of Foundations of Computer Science* 13(1):53–65, 2002. doi:10.1142/S0129054102000947. The DP column of weighted edit distance (search mode, one indel cost) as a solvable automaton, with a vector algorithm that adds once per state but the last; the automaton bitdp's derivation starts from.
- A. Bergeron, S. Hamel. From cascade decompositions to bit-vector algorithms. *Theoretical Computer Science* 313(1):3–16, 2004. doi:10.1016/j.tcs.2003.10.010. Cascade decomposition of solvable automata into threshold levels, each evaluated by one addition; bitdp's Theorem 1 is this cascade for the linear-gap column in closed form.
- G. Navarro. Approximate regular expression searching with arbitrary integer weights. *Nordic Journal of Computing* 11(4):356–373, 2004. Bit-parallel search with integer weights, using counters per automaton state and table lookups.
- L. Allison, T. I. Dix. A bit-string longest-common-subsequence algorithm. *Information Processing Letters* 23(5):305–310, 1986. doi:10.1016/0020-0190(86)90091-8; M. Crochemore, C. S. Iliopoulos, Y. J. Pinzon, J. F. Reid. A fast and practical bit-vector algorithm for the longest common subsequence problem. *Information Processing Letters* 80(6):279–285, 2001. doi:10.1016/S0020-0190(01)00182-X; H. Hyyrö. Bit-parallel LCS-length computation revisited. *AWOCA* 2004, pp. 16–27. Bit-vector LCS, and with it indel distance.
- W. J. Masek, M. S. Paterson. A faster algorithm computing string edit distances. *JCSS* 20(1):18–31, 1980. doi:10.1016/0022-0000(80)90002-1; E. Ukkonen. Finding approximate patterns in strings. *Journal of Algorithms* 6(1):132–137, 1985. doi:10.1016/0196-6774(85)90023-9; S. Wu, U. Manber, G. Myers. A subquadratic algorithm for approximate limited expression matching. *Algorithmica* 15(1):50–67, 1996. doi:10.1007/BF01942606. The DP column read as a finite-state machine.
- T. Mytkowicz, M. Musuvathi, W. Schulte. Data-parallel finite-state machines. *ASPLOS* 2014. doi:10.1145/2541940.2541988. Parallel evaluation of finite-state machines by composing their transition functions, the idea behind reading a DP column as a machine.
- H. Hyyrö, K. Fredriksson, G. Navarro. Increased bit-parallelism for approximate and multiple string matching. *ACM Journal of Experimental Algorithmics* 10, 2005. doi:10.1145/1064546.1180617. Several short patterns packed into one computer word for bit-parallel search, the idea behind `Packed`.

- K. E. Batcher. Sorting networks and their applications. *AFIPS Spring Joint Computer Conference* 32:307–314, 1968. doi:10.1145/1468075.1468121. The odd-even merging network used for unary addition.

## Books

- G. Gopalakrishnan. *Computation Engineering: Applied Automata Theory and Logic*. Springer, 2006. doi:10.1007/0-387-32520-4. Finite-state machines, transducers and boolean function representation.
- D. E. Knuth. *The Art of Computer Programming, Volume 4, Fascicle 6: Satisfiability*. Addison-Wesley, 2015. Cardinality constraints encoded with sorting networks, the source of the idea that merging sorted bit strings counts in unary.
