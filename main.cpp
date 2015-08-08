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

  // cout << "* " << u.rot_max << endl;

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

int put_unit(vector<vector<int>> &bd, const pt &pos, const unit &u, int rot)
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
    for (int x = 0; x < w/2; x++)
      bd[cy][x*2+cy%2] = 0;
    cy--;
  }

  return u.members.size() + 100 * (1 + lines) * lines / 2;
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
  int px = (prob.width - (maxx - minx + 1)) / 2;

  return pt(px*2+(py%2+2)%2, py);
}

double calc_score(const vector<vector<int>> &bd)
{
  int w = bd[0].size(), h = bd.size();

  double height = h, holes = 0, roofs = 0;

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w/2; x++) {
      int xx = x*2 + y%2;
      if (bd[y][xx] != 0) {
        height = min(height, (double)y);
      }

      const int vect[][2] = {
        {-1, -1},
        {1, -1},
        {2, 0},
        {1, 1},
        {-1, 1},
        {-2, 0},
      };

      int rcnt = 0, tcnt = 0;
      for (int i = 0; i < 6; i++) {
        int cx = xx + vect[i][0];
        int cy = y + vect[i][1];

        if (cx >= 0 && cx < w && cy >= 0 && cy < h && bd[cy][cx] == 0) continue;
        if (i < 2) rcnt++;
        tcnt++;
      }

      // if (rcnt == 2) roofs++;
      roofs += rcnt;
      if (tcnt == 6) holes++;
    }
  }

  double ret = 0;

  // vlog() << height << endl;

  ret += height / h * 10000;
  ret += holes * -50;
  ret += roofs * -10;

  return ret;
}

string to_ans(const vector<int> &v)
{
  string ret;
  for (auto &m: v) {
    char c;
    switch (m) {
    case W:   c = 'p'; break;
    case E:   c = 'b'; break;
    case SW:  c = 'a'; break;
    case SE:  c = 'l'; break;
    case CW:  c = 'd'; break;
    case CCW: c = 'k'; break;
    }
    ret += c;
  }
  return ret;
}

struct candidate {
  vector<int> moves;
  int move_score;
  vector<vector<int>> board;
};

candidate gen_cand(vector<vector<int>> bd, const unit &u, pt pos)
{
  vector<int> moves;
  double move_score = 0;
  int rot = 0;

  set<pair<pair<int,int>, int>> ss;

  for (;;) {
    int d = rand() % 6;

    pt prev = pos;
    int prot = rot;

    switch(d) {
    case W:   pos += pt(-2, 0); break;
    case E:   pos += pt(2, 0); break;
    case SW:  pos += pt(-1, 1); break;
    case SE:  pos += pt(1, 1); break;
    case CW:  rot = (rot + 1) % u.rot_max; break;
    case CCW: rot = (rot + u.rot_max - 1) % u.rot_max; break;
    }

    bool ok = check(bd, pos, u, rot);
    if (!ok && (d==CW||d==CCW)) {
      pos = prev;
      rot = prot;
      continue;
    }

    auto key = make_pair(pt2pair(pos), rot);
    if (ss.count(key)) {
      pos = prev;
      rot = prot;
      continue;
    }
    ss.insert(key);

    moves.push_back(d);

    if (!ok) {
      pos = prev;
      rot = prot;
      move_score += put_unit(bd, pos, u, rot);
      break;
    }
  }

  candidate ret;
  ret.moves = moves;
  ret.move_score = move_score;
  ret.board = bd;
  return ret;
}

string solve(const problem &prob, int seed, int tle, int mle)
{
  vector<vector<int>> bd(prob.height, vector<int>(prob.width*2, 1));

  for (int y = 0; y < prob.height; y++)
    for (int x = 0; x < prob.width; x++)
      bd[y][x*2+y%2] = 0;

  for (auto &p: prob.filled)
    bd[p.imag()][p.real()] = 1;

  int move_score = 0;
  vector<int> moves;

  rng r(seed);
  for (int turn = 0; turn < prob.source_length; turn++) {
    vlog() << "turn: " << turn << endl;
    const unit &u = prob.units[r.get() % prob.units.size()];
    pt pos = get_init_pos(prob, u);
    if (!check(bd, pos, u, 0)) break;

    candidate best;
    double best_score = -1e10;
    for (int i = 0; i < 1000; i++) {
      auto cand = gen_cand(bd, u, pos);
      double score = calc_score(cand.board);
      // vlog() << "cand; " << score << ", " << to_ans(cand.first) << endl;
      // print_board(cand.second);
      if (score > best_score) {
        best_score = score;
        best = cand;
      }
    }

    bd = best.board;
    for (auto &m: best.moves)
      moves.push_back(m);
    move_score += best.move_score;

    vlog() << "selected:" << endl;
    print_board(bd);
  }
  cerr << "id: " << prob.id << ", seed: " << seed << ", score: " << move_score << endl;
  print_board(bd);

  return to_ans(moves);
}

int main(int argc, char *argv[])
{
  string tag;
  vector<string> files;
  vector<string> ppw;
  int tle = -1;
  int mle = -1;
  int core = 1;

  for (int i = 1; i < argc; i += 2) {
    string arg = argv[i];

    if (arg == "-f") {
      files.push_back(argv[i+1]);
    }
    if (arg == "-t") {
      tle = atoi(argv[i+1]);
    }
    if (arg == "-m") {
      mle = atoi(argv[i+1]);
    }
    if (arg == "-c") {
      core = atoi(argv[i+1]);
    }
    if (arg == "-p") {
      ppw.push_back(argv[i+1]);
    }
    if (arg == "-g") {
      tag = argv[i+1];
    }

    if (arg == "-v") {
      verbose = true;
    }
  }

  vector<problem> problems;
  for (auto &file: files)
    problems.push_back(read_problem(file));

  int total_problems = 0;
  for (auto &p: problems)
    total_problems += p.source_seeds.size();

  picojson::array v;
  for (auto &p: problems) {
    for (auto &seed: p.source_seeds) {
      string move = solve(p, seed, tle, mle);
      object o;
      o["problemId"] = value((int64_t)p.id);
      o["seed"] = value((int64_t)seed);
      if (tag != "")
        o["tag"] = value(tag);
      o["solution"] = value(move);
      v.push_back(value(o));

      // break;
    }
  }
  cout << value(v) << endl;

  return 0;
}
