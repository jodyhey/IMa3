/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#ifndef _UPDATE_GTREE_H_
#define _UPDATE_GTREE_H_
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS

extern struct genealogy_weights holdgweight_updategenealogy;
extern struct genealogy_weights holdallgweight_updategenealogy;
extern struct probcalc holdallpcalc_updategenealogy;

/* prototype of functions local to update_gtree.cpp */
double findjointime (int ci, int slidepop, int sispop, double edgeuptime,
                     double sisuptime);

double addmigration (int ci, int li);

__END_DECLS
#endif /* _UPDATE_GTREE_H_ */
