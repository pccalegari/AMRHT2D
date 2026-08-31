class dominio {
 private:
  double xbegin, ybegin, xend, yend;
 public:
  dominio ();
  dominio (double xb, double yb, double xe, double ye);
  void set_xbegin(double);
  void set_ybegin(double);
  void set_xend(double);
  void set_yend(double);
  double get_xbegin();
  double get_ybegin();
  double get_xend();
  double get_yend();
};
