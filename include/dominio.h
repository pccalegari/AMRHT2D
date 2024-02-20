class dominio {
 private:
  double xbegin, ybegin, xend, yend;
  int tbc_right, tbc_left, tbc_up, tbc_down;
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
  void set_tbc_left(int tbc);
  void set_tbc_right(int tbc);
  void set_tbc_up(int tbc);
  void set_tbc_down(int tbc);
  int get_tbc_right();
  int get_tbc_left();
  int get_tbc_up();
  int get_tbc_down();
};
