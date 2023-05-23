// https://discourse.mc-stan.org/t/is-it-possible-to-access-the-iteration-step-number-inside-a-stan-program/1871/23
static int itct = 1;

inline void reset_iter(std::ostream* pstream__) {
  itct = 1;
}

inline void add_iter(std::ostream* pstream__) {
  itct += 1;
}
inline int get_iter(std::ostream* pstream__) {
  return itct;
}