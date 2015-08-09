#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <cassert>
#include <set>
#include <unordered_set>

#define PICOJSON_USE_INT64
#include "picojson.h"

using namespace std;
using namespace picojson;

// Mersenne Twister

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}

// main

set<string> ppw;

bool verbose = false;
ostringstream cnull;

ostream &vlog(){
  if (verbose) return cerr;
  return cnull;
}

typedef complex<int> pt;

inline pair<int, int> pt2pair(const pt &p) {
  return make_pair(p.real(), p.imag());
}

inline pt rot_cw(const pt &p)
{
  int x = p.real();
  int y = p.imag();
  return pt((x-3*y)/2, (x + y)/2);
}

inline pt rot_n(const pt &p, int n)
{
  if (n == 0) return p;
  return rot_n(rot_cw(p), n-1);
}

inline pt rot_pivot(const pt &p, const pt &piv, int n)
{
  return rot_n(p - piv, n) + piv;
}

struct unit {
  vector<pt> members;
  pt pivot;

  int rot_max;
  int miny;
};

struct problem {
  int id;
  vector<unit> units;
  int width;
  int height;
  vector<pt> filled;
  int source_length;
  vector<int> source_seeds;
};

ostream &operator<<(ostream &os, const unit &u)
{
  os << "{ " << u.pivot << ", [";
  for (auto &i: u.members)
    os << i << ", ";
  os << "] }";
  return os;
}

ostream &operator<<(ostream &os, const problem &p)
{
  os << "id: " << p.id << endl;
  os << "units:" << endl;
  for (auto &u: p.units)
    os << "  " << u << endl;
  os << "width: " << p.width << endl;
  os << "height: " << p.height << endl;
  os << "filled: ";
  for (auto &q: p.filled)
    os << q << ", ";
  os << endl;
  os << "sourceLength: " << p.source_length << endl;
  os << "sourceSeeds: ";
  for (auto &s: p.source_seeds)
    os << s << ", ";
  os << endl;
  return os;
}

pt to_cell(value &v) {
  int x = v.get<object>()["x"].get<int64_t>();
  int y = v.get<object>()["y"].get<int64_t>();
  return pt(x*2+y%2, y);
}

unit to_unit(value &v) {
  unit u;
  for (auto &c: v.get<object>()["members"].get<picojson::array>()) {
    u.members.push_back(to_cell(c));
  }
  u.pivot = to_cell(v.get<object>()["pivot"]);

  // calculate rotation max
  set<set<pair<int, int>>> ss;
  u.rot_max = 1;
  for (int i = 0; i < 6; i++) {
    set<pair<int, int>> s;
    for (auto &c: u.members) {
      auto p = rot_pivot(c, u.pivot, i);
      s.insert(make_pair(p.real(), p.imag()));
    }
    if (ss.count(s)) break;
    ss.insert(s);
    u.rot_max = i + 1;
  }

  u.miny = 0x7fffffff;
  for (auto &p: u.members)
    u.miny = min(u.miny, p.imag());

  return u;
}

problem to_problem(value &v) {
  problem ret;

  object obj = v.get<object>();
  ret.id = obj["id"].get<int64_t>();
  for (auto &i: obj["units"].get<picojson::array>())
    ret.units.push_back(to_unit(i));
  ret.width = obj["width"].get<int64_t>();
  ret.height = obj["height"].get<int64_t>();
  for (auto &i: obj["filled"].get<picojson::array>())
    ret.filled.push_back(to_cell(i));
  ret.source_length = obj["sourceLength"].get<int64_t>();
  for (auto &i: obj["sourceSeeds"].get<picojson::array>())
    ret.source_seeds.push_back(i.get<int64_t>());

  return ret;
}

problem read_problem(const string &file)
{
  ifstream ifs(file.c_str());
  value v;
  parse(v, ifs);
  return to_problem(v);
}

class rng {
public:
  rng(int seed): a(seed) {}

  int get() {
    int ret = (a >> 16) & 0x7fff;
    a = a * 1103515245 + 12345;
    return ret;
  }

  int a;
};

enum command {
  W, E, SW, SE, CW, CCW
};

char to_char(command m) {
  char c;
  switch (m) {
  // case W:   c = 'p'; break;
  // case E:   c = 'b'; break;
  // case SW:  c = 'a'; break;
  // case SE:  c = 'l'; break;
  // case CW:  c = 'd'; break;
  // case CCW: c = 'k'; break;

  case W:   c = '!'; break;
  case E:   c = 'e'; break;
  case SW:  c = 'i'; break;
  case SE:  c = 'l'; break;
  case CW:  c = 'd'; break;
  case CCW: c = 'k'; break;
  }
  return c;
}

string show_command(command m) {
  switch (m) {
  case W:   return "W";
  case E:   return "E";
  case SW:  return "SW";
  case SE:  return "SE";
  case CW:  return "CW";
  case CCW: return "CCW";
  }
}

map<char, command> cmdmap;

command to_command(char c) {
  c = tolower(c);
  if (cmdmap.empty()) {
    cmdmap['p'] = W;
    cmdmap['\''] = W;
    cmdmap['!'] = W;
    cmdmap['.'] = W;
    cmdmap['0'] = W;
    cmdmap['3'] = W;

    cmdmap['b'] = E;
    cmdmap['c'] = E;
    cmdmap['e'] = E;
    cmdmap['f'] = E;
    cmdmap['y'] = E;
    cmdmap['2'] = E;

    cmdmap['a'] = SW;
    cmdmap['g'] = SW;
    cmdmap['h'] = SW;
    cmdmap['i'] = SW;
    cmdmap['j'] = SW;
    cmdmap['4'] = SW;

    cmdmap['l'] = SE;
    cmdmap['m'] = SE;
    cmdmap['n'] = SE;
    cmdmap['o'] = SE;
    cmdmap[' '] = SE;
    cmdmap['5'] = SE;

    cmdmap['d'] = CW;
    cmdmap['q'] = CW;
    cmdmap['r'] = CW;
    cmdmap['v'] = CW;
    cmdmap['z'] = CW;
    cmdmap['1'] = CW;

    cmdmap['k'] = CCW;
    cmdmap['s'] = CCW;
    cmdmap['t'] = CCW;
    cmdmap['u'] = CCW;
    cmdmap['w'] = CCW;
    cmdmap['x'] = CCW;
  }

  return cmdmap[c];
}

inline bool check_one(const vector<vector<int>> &bd,
                      const pt &pos, const pt &pivot, const pt &cell, int rot)
{
  pt e = rot_pivot(cell + pos, pivot + pos, rot);
  int w = bd[0].size(), h = bd.size();
  return e.real() >= 0 && e.real() < w && e.imag() >= 0 && e.imag() < h && bd[e.imag()][e.real()] == 0;
}

bool check(const vector<vector<int>> &bd, const pt &pos, const unit &u, int rot)
{
  for (auto &c: u.members)
    if (!check_one(bd, pos, u.pivot, c, rot))
      return false;
  return true;
}

pair<int,int> put_unit(vector<vector<int>> &bd, const pt &pos, const unit &u, int rot, int ls_old)
{
  for (auto &c: u.members) {
    pt p = pos + rot_pivot(c, u.pivot, rot);
    assert(p.imag()>=0&&p.imag()<bd.size()&&p.real()>=0&&p.real()<bd[0].size());
    bd[p.imag()][p.real()] = 1;
  }

  int w = bd[0].size()/2, h = bd.size();
  int cy = h - 1;
  int lines = 0;

  for (int y = h - 1; y >= 0; y--) {
    bool all = true;
    for (int x = 0; x < w; x++) {
      if (bd[y][x*2+y%2] == 0) {
        all = false;
        break;
      }
    }

    if (all) {
      lines++;
    }
    else {
      if (cy != y) {
        for (int x = 0; x < w; x++)
          bd[cy][x*2+cy%2] = bd[y][x*2+y%2];
      }
      cy--;
    }
  }

  while (cy >= 0) {
    for (int x = 0; x < w; x++)
      bd[cy][x*2+cy%2] = 0;
    cy--;
  }

  int points = u.members.size() + 100 * (1 + lines) * lines / 2;
  int bonus = ls_old > 1 ? points * (ls_old - 1) / 10 : 0;

  return make_pair(points + bonus, lines);
}

void print_board(const vector<vector<int>> &bd)
{
  int w = bd[0].size()/2, h = bd.size();
  for (int y = 0; y < h; y++) {
    for (int xx = 0; xx < w; xx++) {
      int x = xx * 2 + y % 2;
      if (y % 2 == 1) vlog() << " ";
      vlog() << (bd[y][x]==1 ? "o" : bd[y][x]==2 ? "#" : bd[y][x]==3 ? "@" : ".");
      if (y % 2 == 0) vlog() << " ";
    }
    vlog() << endl;
  }
  vlog() << endl;
}

void print_board_(const vector<vector<int>> &bd_, const pt &pos, const unit &u, int rot)
{
  vector<vector<int>> bd = bd_;

  for (auto &c: u.members) {
    pt p = pos + rot_pivot(c, u.pivot, rot);
    assert(p.imag()>=0&&p.imag()<bd.size()&&p.real()>=0&&p.real()<bd[0].size());
    bd[p.imag()][p.real()] = 2;
  }
  pt q = pos + u.pivot;
  if (q.real() >= 0 && q.real() < (int)bd[0].size() &&
      q.imag() >= 0 && q.imag() < (int)bd.size())
    bd[q.imag()][q.real()] = 3;
  print_board(bd);
}

pt get_init_pos(const problem &prob, const unit &u)
{
  int maxx = -0x7fffffff, minx = 0x7fffffff, miny = 0x7fffffff;

  for (auto &p: u.members) {
    maxx = max(maxx, p.real()/2);
    minx = min(minx, p.real()/2);
    miny = min(miny, p.imag());
  }

  int py = -miny;
  int px = (prob.width - (maxx - minx + 1)) / 2 - minx;

  return pt(px*2+(py%2+2)%2, py);
}

vector<double> eparam;

const int FN = 9;
const int NUM_FEATURES = FN+3;

double calc_score(const vector<vector<int>> &bd, const pt &last, int lines)
{
  int w = bd[0].size(), h = bd.size();

  double height = h; //, holes = 0, roofs = 0, hroofs = 0, sides = 0, hsides = 0;

  double feat[FN] = {};

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w/2; x++) {
      int xx = x*2 + y%2;
      if (bd[y][xx] != 0) {
        height = min(height, (double)y);
      }

      if (bd[y][xx] != 0)
        continue;

      const int vect[][2] = {
        {-1, -1},
        {1, -1},
        {2, 0},
        {-2, 0},
        {1, 1},
        {-1, 1},
      };

      int pat = 0;
      for (int i = 0; i < 6; i++) {
        int cx = xx + vect[i][0];
        int cy = y + vect[i][1];

        if (!(cx >= 0 && cx < w && cy >= 0 && cy < h && bd[cy][cx] == 0)) continue;

        pat |= (1 << i);
      }

      int f = 0;
      f = f * 3 + ((pat>>0)&1) + ((pat>>1)&1);
      f = f * 3 + ((pat>>2)&1) + ((pat>>3)&1);
      // f = f * 3 + ((pat>>4)&1) + ((pat>>5)&1);
      // f = f * 2 + (bd[y][xx] != 0);

      feat[f]++;
    }
  }

  double ret = 0;

  // vlog() << height << ", " << holes << ", " << roofs << ", " << sides << endl;
  // height /= h;
  // if (height >= 0.8) height = 1;

  for (int i = 0; i < FN; i++)
    ret += eparam[i] * feat[i];

  ret += eparam[FN + 0] * height * w;
  ret += eparam[FN + 1] * (double)last.imag() * w;
  ret += eparam[FN + 2] * (double)(lines+1)*lines/2 * w * h;

  return ret;
}

string to_ans(const vector<int> &v)
{
  string ret;
  for (auto &m: v)
    ret += to_char((command)m);
  return ret;
}

struct candidate {
  candidate() {
    score = -1e10;
  }

  candidate(const candidate &r)
    : moves(r.moves)
    , move_score(r.move_score)
    , ls_old(r.ls_old)
    , board(r.board)
    , score(r.score)
  { }

  candidate &operator=(const candidate &r) {
    this->moves = r.moves;
    this->move_score = r.move_score;
    this->board = r.board;
    this->score = r.score;
    this->ls_old = r.ls_old;
    return *this;
  }

  vector<int> moves;
  int move_score;
  int ls_old;
  vector<vector<int>> board;
  double score;

  bool operator<(const candidate &r) const {
    return score < r.score;
  }
};

// candidate gen_cand(vector<vector<int>> bd, const unit &u, pt pos)
// {
//   vector<int> moves;
//   double move_score = 0;
//   int rot = 0;

//   set<pair<pair<int,int>, int>> ss;

//   for (;;) {
//     int d = rand() % 6;

//     pt prev = pos;
//     int prot = rot;

//     switch(d) {
//     case W:   pos += pt(-2, 0); break;
//     case E:   pos += pt(2, 0); break;
//     case SW:  pos += pt(-1, 1); break;
//     case SE:  pos += pt(1, 1); break;
//     case CW:  rot = (rot + 1) % u.rot_max; break;
//     case CCW: rot = (rot + u.rot_max - 1) % u.rot_max; break;
//     }

//     bool ok = check(bd, pos, u, rot);
//     if (!ok && (d==CW||d==CCW)) {
//       pos = prev;
//       rot = prot;
//       continue;
//     }

//     auto key = make_pair(pt2pair(pos), rot);
//     if (ss.count(key)) {
//       pos = prev;
//       rot = prot;
//       continue;
//     }
//     ss.insert(key);

//     moves.push_back(d);

//     if (!ok) {
//       pos = prev;
//       rot = prot;
//       move_score += put_unit(bd, pos, u, rot);
//       break;
//     }
//   }

//   candidate ret;
//   ret.moves = moves;
//   ret.move_score = move_score;
//   ret.board = bd;
//   return ret;
// }

void rec(const vector<vector<int>> &bd, const unit &u, const pt &pos, int rot,
         vector<int> &hist, set<pair<pair<int,int>, int>> &ss, candidate &cand, int ls_old, int prev)
{
  auto key = make_pair(make_pair(pos.real(), pos.imag()), rot);
  if (ss.count(key)) return;
  ss.insert(key);

  int invalid = -1;

  static const int tbl[][6] = {
    {1,0,2,3,4,5},
    {2,1,0,3,4,5},
    {0,1,2,3,4,5}
  };

  const int *ord = tbl[2];

  // const int *ord = tbl[0];
  // if (prev == E) ord = tbl[1];
  // if (prev == SW) ord = tbl[2];

  // int ord[] = {0, 1, 2, 3, 4, 5};
  // for (int i = 0; i < 6; i++)
  //   swap(ord[i], ord[i+genrand_int31()%(6-i)]);

  for (int mov_ = 0; mov_ < 6; mov_++) {
    int mov = ord[mov_];
    pt npos = pos;
    int nrot = rot;
    switch(mov) {
    case W:   npos += pt(-2, 0); break;
    case E:   npos += pt(2, 0); break;
    case SW:  npos += pt(-1, 1); break;
    case SE:  npos += pt(1, 1); break;
    case CW:  nrot = (nrot + 1) % u.rot_max; break;
    case CCW: nrot = (nrot + u.rot_max - 1) % u.rot_max; break;
    }
    if (!check(bd, npos, u, nrot)) {
      invalid = mov;
      continue;
    }

    hist.push_back(mov);
    rec(bd, u, npos, nrot, hist, ss, cand, ls_old, mov);
    hist.pop_back();
  }

  if (invalid >= 0) {
    candidate tmp;
    hist.push_back(invalid);
    tmp.moves = hist;
    hist.pop_back();

    vector<vector<int>> bdd = bd;
    auto tt = put_unit(bdd, pos, u, rot, ls_old);
    tmp.move_score = tt.first;
    tmp.board = bdd;

    pt upper = pt(0, 0x7fffffff);
    for (auto &p: u.members) {
      pt q = rot_pivot(p + pos, u.pivot + pos, rot);
      if (q.imag() < upper.imag())
        upper = q;
    }

    tmp.ls_old = tt.second;

    tmp.score = calc_score(bdd, upper, tt.second);
    cand = max(cand, tmp);
  }

  return;
}

vector<vector<int>> make_board(const problem &prob)
{
  vector<vector<int>> bd(prob.height, vector<int>(prob.width*2, 1));

  for (int y = 0; y < prob.height; y++)
    for (int x = 0; x < prob.width; x++)
      bd[y][x*2+y%2] = 0;

  for (auto &p: prob.filled)
    bd[p.imag()][p.real()] = 1;

  return bd;
}

int power_score(const string &cmds)
{
  // return 0;

  int ret = 0;

  for (auto &word: ppw) {
    int reps = 0;
    for (int i = 0; i < (int)cmds.size(); i++) {
      bool ok = true;
      for (int j = 0; j < (int)word.size(); j++) {
        if (i + j < (int)cmds.size() && tolower(cmds[i+j]) != tolower(word[j])) {
          ok = false;
          break;
        }
      }
      if (ok) reps++;
    }

    if (reps > 0) ret += 300;
    ret += 2 * word.length() * reps;
  }

  return ret;
}

pair<string, int> solve(const problem &prob, int seed, int tle, int mle, const vector<double> &param, bool progress)
{
  eparam = param;

  auto bd = make_board(prob);
  if (progress) print_board(bd);

  int move_score = 0;
  vector<int> moves;

  int ls_old = 0;

  rng r(seed);
  for (int turn = 0; turn < prob.source_length; turn++) {
    // vlog() << "turn: " << turn << endl;
    const unit &u = prob.units[r.get() % prob.units.size()];
    pt pos = get_init_pos(prob, u);
    if (!check(bd, pos, u, 0)) break;

    vector<int> hist;
    set<pair<pair<int,int>, int>> ss;
    candidate best;
    rec(bd, u, pos, 0, hist, ss, best, ls_old, CW);

    bd = best.board;
    for (auto &m: best.moves)
      moves.push_back(m);
    move_score += best.move_score;
    ls_old = best.ls_old;

    if (progress) {
      vlog() << "selected: " << best.score << endl;
      print_board(bd);
    }
  }
  if (progress) {
    cerr << "id: " << prob.id << ", seed: " << seed << ", score: " << move_score << endl;
    print_board(bd);
  }

  return make_pair(to_ans(moves), move_score + power_score(to_ans(moves)));
}

void output_solution(int problem_id, int seed, int score, const string &moves, const string &tag)
{
  stringstream ss;
  ss << "out-full/" << problem_id << "-" << seed << "-" << score << ".json";
  ofstream ofs(ss.str().c_str());

  object o;
  o["problemId"] = value((int64_t)problem_id);
  o["seed"] = value((int64_t)seed);
  if (tag != "")
    o["tag"] = value(tag);
  o["solution"] = value(moves);

  ofs << value(o) << endl;
}

pair<string, int> annealing(const problem &p, int seed, int tle, int mle, const string &tag,
                            double init_temp, double temp_decay, int rounds)
{
  int best_score = -1;
  string best_move;

  vector<double> pp(NUM_FEATURES);
  for (int i = 0; i < (int)pp.size(); i++)
    pp[i] = genrand_real1() - 0.5;
  double score = -1;

  double temp = init_temp;
  for (int t = 0; t < rounds; t++, temp *= temp_decay) {
    if (t % 10 == 0)
      cerr << "round: [" << t << "/" << rounds << "]: "
           << temp << ", " << score << "        \r" << flush;

    if (temp < 1){
      cerr << "REBOOT                                   " << endl;
      temp = best_score / 10;
    }

    auto bkup = pp;

    if (genrand_int31() % 2 == 0) {
      pp[genrand_int31() % pp.size()] = genrand_real1() - 0.5;
    }
    else {
      // pp[genrand_int31() % pp.size()] *= genrand_real1() * 0.5 - 0.25 + 1.0;
      pp[genrand_int31() % pp.size()] = genrand_real1() - 0.5;
    }

    auto rr = solve(p, seed, tle, mle, pp, false);
    if (score <= rr.second ||
        genrand_real2() <= exp((rr.second - score) / temp)) {
      score = rr.second;

      if (best_score < rr.second) {
        best_score = rr.second;
        best_move = rr.first;

        cerr << "*BEST* " << p.id << "-" << seed << ", round " << t << ": " << best_score
             << "                      " << endl;
        // cerr << "param[] = [";
        // for (int i = 0; i < pp.size(); i++)
        //   cerr << pp[i] << ", ";
        // cerr << "]" << endl;;

        output_solution(p.id, seed, best_score, best_move, tag);
      }
    }
    else {
      pp = bkup;
    }
  }

  return make_pair(best_move, best_score);
}

vector<double> select(vector<pair<pair<int, string>, vector<double>>> &ord)
{
  double total_adapt = 0;
  for (auto &c: ord)
    total_adapt += pow(log((double)c.first.first), 2);

  for (auto &c: ord) {
    double adp = pow(log((double)c.first.first), 2);
    if (genrand_real1() < (adp/total_adapt))
      return c.second;
    total_adapt -= adp;
  }

  return ord[genrand_int31() % ord.size()].second;
}

pair<string, int> ga(const problem &p, int seed, int tle, int mle, const string &tag)
{
  int best_score = -1;
  string best_move;

  int num_cands = 500;
  vector<vector<double>> cands;

  for (int i = 0; i < num_cands; i++) {
    vector<double> pp(NUM_FEATURES);
    for (int i = 0; i < (int)pp.size(); i++)
      pp[i] = genrand_real1();
    cands.push_back(pp);
  }

  for (int gen = 0; gen < 10000; gen++) {
    vector<pair<pair<int, string>, vector<double>>> ord;
    for (auto &ncand: cands) {
      auto rr = solve(p, seed, tle, mle, ncand, false);
      ord.emplace_back(make_pair(rr.second, rr.first), ncand);
    }

    int cur_best = -1, cur_best_ix = 0;
    for (int i = 0; i < num_cands; i++) {
      if (ord[i].first.first > cur_best) {
        cur_best = ord[i].first.first;
        cur_best_ix = i;
      }
    }
    cerr << "GEN: " << gen << ", score = " << cur_best << endl;

    if (best_score < cur_best) {
      best_score = ord[cur_best_ix].first.first;
      best_move = ord[cur_best_ix].first.second;

      cerr << "*BEST* " << p.id << "-" << seed << ": " << best_score << endl;
      output_solution(p.id, seed, best_score, best_move, tag);
    }

    vector<vector<double>> next_cands;

    while ((int)next_cands.size() < num_cands) {
      double x = genrand_real1();

      if (x < 0.1) {
        next_cands.push_back(select(ord));
      }
      else if (x < 0.99) {
        if (num_cands - next_cands.size() < 2) continue;

        vector<double> sa = select(ord);
        vector<double> sb = select(ord);

        for (int i = 0; i < (int)sa.size(); i++) {
          if (genrand_real1() < 0.5)
            swap(sa[i], sb[i]);
        }
        next_cands.push_back(sa);
        next_cands.push_back(sb);
      }
      else {
        vector<double> sa = select(ord);
        sa[genrand_int31()%sa.size()] = genrand_real1();
      }
    }

    cands = next_cands;
  }

  return make_pair(best_move, best_score);
}

void replay(const problem &prob, int seed, const string &moves)
{
  vlog() << "replay:" << endl;
  auto bd = make_board(prob);
  rng r(seed);

  const char *p = moves.c_str();
  int move_score = 0;

  for (int turn = 0; turn < prob.source_length; turn++) {
    const unit &u = prob.units[r.get() % prob.units.size()];
    pt pos = get_init_pos(prob, u);
    int rot = 0;
    if (!check(bd, pos, u, rot))
      break;
    vlog() << "*POP*" << endl;
    print_board_(bd, pos, u, rot);

    for (;*p;) {
      pt npos = pos;
      int nrot = rot;
      command mov = to_command(*p++);
      switch(mov) {
      case W:   npos += pt(-2, 0); break;
      case E:   npos += pt(2, 0); break;
      case SW:  npos += pt(-1, 1); break;
      case SE:  npos += pt(1, 1); break;
      case CW:  nrot = (nrot + 1) % u.rot_max; break;
      case CCW: nrot = (nrot + u.rot_max - 1) % u.rot_max; break;
      }
      if (check(bd, npos, u, nrot)) {
        pos = npos;
        rot = nrot;
        vlog() << "cmd: " << show_command(mov) << endl;
        print_board_(bd, pos, u, rot);
      }
      else {
        vlog() << "cmd: LOCK (" << show_command(mov) << ")" << endl;
        move_score += put_unit(bd, pos, u, rot, 0).first;
        print_board(bd);
        break;
      }
    }
  }

  auto psc = power_score(moves);

  vlog() << "move score: " << move_score << endl;
  vlog() << "power score: " << psc << endl;
  vlog() << "total score: " << move_score + psc << endl;
}

int main(int argc, char *argv[])
{
  ppw.insert("ei!");
  ppw.insert("yuggoth");

  // srand(time(NULL));
  init_genrand(time(NULL));

  string tag;
  vector<string> files;
  int tle = -1;
  int mle = -1;
  int core = 1;

  bool anneal = false;
  double init_temp;
  double temp_decay;
  int turns;

  bool gene = false;

  vector<double> def_param {-2.016, 2.928, 4.761, 4.257, 1.862, -0.857, 1.512,};

  for (int i = 1; i < argc; ) {
    string arg = argv[i];

    if (arg == "-f") {
      files.push_back(argv[i+1]);
      i += 2;
    }
    else if (arg == "-t") {
      tle = atoi(argv[i+1]);
      i += 2;
    }
    else if (arg == "-m") {
      mle = atoi(argv[i+1]);
      i += 2;
    }
    else if (arg == "-c") {
      core = atoi(argv[i+1]);
      i += 2;
    }
    else if (arg == "-p") {
      string w = argv[i+1];
      for (auto &c: w) c = tolower(c);
      ppw.insert(w);
      i += 2;
    }

    else if (arg == "-g") {
      tag = argv[i+1];
      i += 2;
    }
    else if (arg == "-v") {
      verbose = true;
      i++;
    }
    else if (arg == "-a") {
      anneal = true;

      istringstream iss(argv[i+1]);
      char dmy;
      iss >> init_temp >> dmy >> temp_decay >> dmy >> turns;

      i += 2;
    }
    else if (arg == "-ga") {
      gene = true;
      i++;
    }
  }

  vector<problem> problems;
  for (auto &file: files)
    problems.push_back(read_problem(file));

  int total_problems = 0;
  for (auto &p: problems)
    total_problems += p.source_seeds.size();

  picojson::array v;
  vector<vector<int>> scores;
  for (auto &p: problems) {
    vector<int> ss;
    cerr << "\n" << "problem: " << p.id << ": " << p.width << "x" << p.height << endl;

    int ttt = 0;
    for (auto &seed: p.source_seeds) {
      cerr << "\n" << p.id << " [" << ++ttt << "/" << p.source_seeds.size() << "]" << endl;

      auto sol =
        anneal ? annealing(p, seed, tle, mle, tag, init_temp, temp_decay, turns):
        gene ? ga(p, seed, tle, mle, tag):
        solve(p, seed, tle, mle, def_param, true);

      if (!anneal) replay(p, seed, sol.first);

      object o;
      o["problemId"] = value((int64_t)p.id);
      o["seed"] = value((int64_t)seed);
      if (tag != "")
        o["tag"] = value(tag);
      o["solution"] = value(sol.first);
      v.push_back(value(o));
      ss.push_back(sol.second);
    }
    scores.push_back(ss);
  }
  cout << value(v) << endl;

  double gtot = 0;

  cerr << "summary: " << endl;
  vector<pair<int,double>> vv;
  for (int i=0;i<(int)problems.size();i++) {
    double avg = 0;
    for (auto &p: scores[i]) avg += p;
    avg /= scores[i].size();
    vv.emplace_back(problems[i].id, avg);

    double avg_size = 0;
    for (auto &p: problems[i].units)
      avg_size += p.members.size();
    avg_size /= problems[i].units.size();
    double blks = problems[i].source_length * avg_size;
    double exp_line = blks / problems[i].width;
    gtot += avg / exp_line;
  }

  sort(vv.begin(), vv.end());

  for (auto &i: vv)
    cerr << i.first << ": " << i.second << endl;

  cerr << "grand total: " << gtot / problems.size() << endl;

  return 0;
}
