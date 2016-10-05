#include <Siv3D.hpp>
#include <stack>
#include <queue>

using namespace std;

class ExtraPiece {
public:

	const Array<Point> getPiece(const Image& image, const double max_distance = 3.0) {

		const Polygon polygon = Imaging::FindContour(image, false);
		const Array<Vec2> vecArray = polygon.simplified(max_distance).outer();

		Array<Point> vec;
		for (size_t i = 0, size = vecArray.size(); i < size; i++)
		{
			const double range = Math::Sqrt(Math::Square(vecArray[i].x - vecArray[(i + 1) % size].x) + Math::Square(vecArray[i].y - vecArray[(i + 1) % size].y));
			const double angel = Math::Atan2(vecArray[(i + 1) % size].y - vecArray[i].y, vecArray[(i + 1) % size].x - vecArray[i].x);

			//最低の辺の長さ未満ならば辺を削除する
			if (range > 10)
			{
				vec.push_back(Point(vecArray[i].x, vecArray[i].y));
			}
			else
			{
				const Vec2 p1 = vecArray[(i + size - 1) % size];
				const Vec2 p2 = vecArray[(i + 0) % size];
				const Vec2 p3 = vecArray[(i + 1) % size];
				const Vec2 p4 = vecArray[(i + 2) % size];

				const Vec2 v1 = p2 - p1;
				const Vec2 v2 = p3 - p4;

				const double sinv1 = (v1.x != 0) ? (v1.y / v1.x) : (1);
				const double sinv2 = (v2.x != 0) ? (v2.y / v2.x) : (1);

				const double x1 = ((p1.x*sinv1 - p4.x*sinv2) - (p1.y - p4.y)) / (sinv1 - sinv2);
				const double y1 = p1.y + sinv1*(x1 - p1.x);
				const double x2 = ((p1.x*sinv1 - p4.x*sinv2) - (p1.y - p4.y)) / (sinv1 - sinv2);
				const double y2 = p4.y + sinv2*(x2 - p4.x);

				vec.push_back(Point((int)x1, (int)y1));
			}
		}

		return vec;
	}

};

class ImageChip {
public:

	static const Array<Image> chip(Image image) {

		Grid<int> bmp(image.size, 0);
		Grid<int> flag(image.size, -1);

		for (const auto& p : step(image.size))
			bmp[p] = image[p].r;

		const int pieceNumber = numbering(bmp, flag);

		return chiper(bmp, flag);
	}

private:

	static const inline bool inside(const Point& p, const Size& s) {
		return (0 <= p.x && 0 <= p.y && p.x < s.x && p.y < s.y);
	}

	static const int numbering(const Grid<int>& bmp, Grid<int>& flag) {

		const auto inside = [](const Point& p, const Size& s) {
			return (0 <= p.x && 0 <= p.y && p.x < s.x && p.y < s.y);
		};

		int pieceCount = 0;
		for (const auto& p : step(bmp.size()))
		{
			if (flag[p] != -1) continue;

			std::stack<Point> sta;

			sta.push(p);
			const int color = bmp[p];
			flag[p] = pieceCount;

			while (!sta.empty())
			{
				const auto point = sta.top();
				sta.pop();

				for (int dy = -1; dy <= 1; dy++)
				{
					for (int dx = -1; dx <= 1; dx++)
					{
						if (dy == 0 && dx == 0) continue;

						const Point dp = point + Point(dx, dy);

						if (inside(dp, bmp.size()) && flag[dp] == -1 && bmp[dp] == color)
						{
							sta.push(dp);
							flag[dp] = pieceCount;
						}

					}
				}

			}
			pieceCount++;
		}
		/*
		Image img(bmp.size(), Palette::Black);
		for (const auto& p : step(bmp.size()))
		{
		const int c = 255 / pieceCount*flag[p];
		img[p] = Color(c, c, c);
		}
		img.savePNG(L"切り分け.png");
		*/
		return pieceCount;
	}

	static const Array<Image> chiper(const Grid<int>& bmp, const Grid<int>& label) {

		Array<Image> images;

		int pieceCount = 1;
		Grid<bool> flag(bmp.size(), false);

		for (const auto& p : step(bmp.size()))
		{
			if (label[p] <= 0) continue;

			if (label[p] == pieceCount && !flag[p])
			{
				Array<Point> points;

				std::stack<Point> sta;
				sta.push(p);
				flag[p] = true;
				points.push_back(p);

				int left = p.x, right = p.x, up = p.y, down = p.y;

				while (!sta.empty())
				{
					const Point point = sta.top();
					sta.pop();

					for (int dy = -1; dy <= 1; dy++)
					{
						for (int dx = -1; dx <= 1; dx++)
						{
							if (dy == 0 && dx == 0) continue;

							const Point dp = point + Point(dx, dy);

							if (inside(dp, bmp.size()) && !flag[dp] && label[dp] == pieceCount)
							{
								sta.push(dp);
								points.push_back(dp);
								flag[dp] = true;

								left = std::min(left, dp.x);
								right = std::max(right, dp.x);
								up = std::min(up, dp.y);
								down = std::max(down, dp.y);
							}

						}
					}

				}

				Image chipImage({ right - left + 1,down - up + 1 }, Palette::Black);

				for (const auto& point : points)
				{
					chipImage[point - Point(left, up)] = Palette::White;
				}
				images.push_back(chipImage);
				pieceCount++;
			}

		}

		return images;
	}

};

class Filter {
public:

	static Image filter(Image image) {

		Grid<bool> check(image.size, false);

		for (const auto& p : step(image.size))
		{
			if (!check[p])
			{
				check[p] = true;

				std::queue<Point> que;
				que.push(p);

				Array<Point> points;
				points.push_back(p);

				const auto inside = [](const Point& p, const Size& s) { return (0 <= p.x && p.x < s.x && 0 <= p.y && p.y < s.y); };

				while (!que.empty())
				{
					const Point point = que.front();
					que.pop();

					const Point p1 = point + Point(0, -1);
					const Point p2 = point + Point(-1, 0);
					const Point p3 = point + Point(1, 0);
					const Point p4 = point + Point(0, 1);

					if (inside(p1, image.size) && !check[p1] && image[p1] == Palette::White)
					{
						check[p1] = true;
						que.push(p1);
						points.push_back(p1);
					}
					if (inside(p2, image.size) && !check[p2] && image[p2] == Palette::White)
					{
						check[p2] = true;
						que.push(p2);
						points.push_back(p2);
					}
					if (inside(p3, image.size) && !check[p3] && image[p3] == Palette::White)
					{
						check[p3] = true;
						que.push(p3);
						points.push_back(p3);
					}
					if (inside(p4, image.size) && !check[p4] && image[p4] == Palette::White)
					{
						check[p4] = true;
						que.push(p4);
						points.push_back(p4);
					}
				}

				if (points.size() < 200)
				{
					for (const auto& pos : points)
					{
						image[pos] = Palette::Black;
					}
				}

			}


		}

		return image;
	}

	static const Image median(Image image) {

		Grid<int> bmp(image.size, 0);
		Grid<int> bmp2(image.size, 0);

		for (const auto& p : step(image.size))
			bmp[p] = image[p].r;

		const int size = 3;

		for (int y = 0; y < image.size.y; y++)
		{
			for (int x = 0; x < image.size.x; x++)
			{
				int val[size*size];
				int count = 0;
				for (int dy = -size / 2; dy <= size / 2; dy++)
				{
					for (int dx = -size / 2; dx <= size / 2; dx++)
					{
						const Point p2 = Point(x + dx, y + dy);
						if (inside(p2, image.size))
						{
							val[count++] = bmp[p2];
						}
					}
				}
				std::sort(val, val + count);
				bmp2[y][x] = val[count / 2];
			}
		}

		for (int y = 0; y < image.size.y; y++)
		{
			for (int x = 0; x < image.size.x; x++)
			{
				const int val = bmp2[y][x];
				image[y][x] = Color(val, val, val);
			}
		}

		return image;
	}

private:
	static const inline bool inside(const Point& p, const Size& s) {
		return (0 <= p.x && 0 <= p.y && p.x < s.x && p.y < s.y);
	}

};

struct Corner {
	Point point;
	double length;
	double angle;
};

struct Piece {
	Array<Corner> corners;
};

class AI {
public:

	void think(const Array<Array<Point>>& points) {

		getData(points);

		const size_t size = pieces.size();

		struct PairMatch {
			PairMatch() = default;
			PairMatch(size_t _n1, size_t _n2, int s) : n1(_n1), n2(_n2), score(s) {}

			size_t n1;
			size_t n2;
			int score;

			const bool operator<(const PairMatch& other) const { return score < other.score; }

		};

		priority_queue<PairMatch> scoreQue;
		vector<PairMatch> scoreVec;

		for (size_t i = 0; i < size; i++)
		{
			for (size_t j = i + 1; j < size; j++)
			{
				const int score = match(i, j);
				if (score > 0)
					scoreVec.push_back(PairMatch{ i, j, score });
				scoreQue.push(PairMatch{ i, j, score });
			}
		}

		for (const auto& ans : scoreVec)
		{
			printf("%3d, %3d: %d\n", (int)ans.n1, (int)ans.n2, ans.score);
		}
		cout << endl;
		while (!scoreQue.empty())
		{
			const auto data = scoreQue.top();
			scoreQue.pop();
			printf("%3d, %3d: %d\n", (int)data.n1, (int)data.n2, data.score);
		}

	}

private:

	Array<Piece> pieces;

	void getData(const Array<Array<Point>>& points) {

		for (const auto& piece : points)
		{
			Array<double> length;
			Array<double> angle;

			Piece pie;
			Corner corner;

			for (size_t i = 0, size = piece.size(); i < size; i++)
			{
				const Point p1 = piece[(i + 0) % size];
				const Point p2 = piece[(i + 1) % size];

				const Point _p1 = piece[(i + size - 1) % size];
				const Point _p2 = piece[(i + size - 0) % size];

				const Point dp = p2 - p1;
				const Point _dp = _p2 - _p1;

				const double len = Math::Sqrt(Math::Square(p1.x - p2.x) + Math::Square(p1.y - p2.y));

				const double ang = Math::Atan2((double)dp.x, (double)dp.y);
				const double _ang = Math::Atan2((double)_dp.x, (double)_dp.y);

				const double dang = Pi - (_ang - ang);

				corner.point = piece[i];
				corner.length = len;
				corner.angle = dang;
				pie.corners.push_back(corner);
			}
			pieces.emplace_back(pie);
		}

	}

	const int match(const size_t n1, const size_t n2) const {

		const double angleErrer = 1.0 / 180 * Pi;
		const double lengthErrer = 3;

		const auto size1 = pieces[n1].corners.size();
		const auto size2 = pieces[n2].corners.size();

		int maxScore = 0;

		for (const auto i : step(size1))
		{
			for (const auto j : step(size2))
			{
				const double addAngle = pieces[n1].corners[i].angle + pieces[n2].corners[j].angle;

				const double len1_1 = pieces[n1].corners[i].length;
				const double len1_2 = pieces[n2].corners[j].length;

				const double len2_1 = pieces[n1].corners[(i + size1 - 1) % size1].length;
				const double len2_2 = pieces[n2].corners[(j + size2 - 1) % size2].length;

				const double ang1_1 = pieces[n1].corners[(i + 1) % pieces[n1].corners.size()].angle;
				const double ang1_2 = pieces[n2].corners[(j + 1) % pieces[n2].corners.size()].angle;

				const double ang2_1 = pieces[n1].corners[(i + size1 - 1) % size1].angle;
				const double ang2_2 = pieces[n2].corners[(j + size2 - 1) % size2].angle;

				int score = 0;

				//角の合計が180度
				if (abs(addAngle - Pi) <= angleErrer)
				{
					//len1が長い場合
					if (len1_1 - len1_2 + lengthErrer > 0)
					{
						//len2の次の角度が180度を越えない
						if (ang1_2 < Pi + angleErrer)
						{
							score += 100;
						}
					}
					//len2が長い場合
					else if (len1_2 - len1_1 + lengthErrer > 0)
					{
						if (ang1_1 < Pi + angleErrer)
						{
							score += 100;
						}
					}
					//辺の長さが同じ
					else
					{
						//角度の合計が180度
						if (abs(ang1_1 + ang1_2 - Pi) < angleErrer)
						{
							score += 200;
						}
						//角度の合計が360度
						if (abs(ang1_1 + ang1_2 - Pi * 2) < angleErrer)
						{
							score += 400;
						}
					}

				}

				//角の合計が360度
				if (abs(addAngle - Pi * 2) <= angleErrer)
				{
					//len1が長い場合
					if (len1_1 - len1_2 + lengthErrer > 0)
					{
						//len2の次の角度が180度を越えない
						if (ang1_2 < Pi + angleErrer)
						{
							if (len2_1 - len2_2 + lengthErrer > 0)
							{
								if (ang2_1 < Pi + angleErrer)
								{
									score += 400;
								}
							}
							else if (len2_2 - len2_1 + lengthErrer > 0)
							{
								if (ang2_2 < Pi + angleErrer)
								{
									score += 400;
								}
							}
							else
							{
								//角度の合計が180度
								if (abs(ang2_1 + ang2_2 - Pi) < angleErrer)
								{
									score += 200;
								}
								//角度の合計が360度
								if (abs(ang2_1 + ang2_2 - Pi * 2) < angleErrer)
								{
									score += 400;
								}
							}
						}
					}
					//len2が長い場合
					else if (len1_2 - len1_1 + lengthErrer > 0)
					{
						//len1の次の角度が180度を越えない
						if (ang1_1 < Pi + angleErrer)
						{
							if (len2_1 - len2_2 + lengthErrer > 0)
							{
								if (ang2_1 < Pi + angleErrer)
								{
									score += 400;
								}
							}
							else if (len2_2 - len2_1 + lengthErrer > 0)
							{
								if (ang2_2 < Pi + angleErrer)
								{
									score += 400;
								}
							}
							else
							{
								//角度の合計が180度
								if (abs(ang2_1 + ang2_2 - Pi) < angleErrer)
								{
									score += 200;
								}
								//角度の合計が360度
								if (abs(ang2_1 + ang2_2 - Pi * 2) < angleErrer)
								{
									score += 400;
								}
							}
						}
					}
					//辺の長さが同じ
					else
					{
						//角度の合計が180度
						if (abs(ang1_1 + ang1_2 - Pi) < angleErrer)
						{
							if (len2_1 - len2_2 + lengthErrer > 0)
							{
								if (ang2_1 < Pi + angleErrer)
								{
									score += 800;
								}
							}
							else if (len2_2 - len2_1 + lengthErrer > 0)
							{
								if (ang2_2 < Pi + angleErrer)
								{
									score += 800;
								}
							}
							else
							{
								//角度の合計が180度
								if (abs(ang2_1 + ang2_2 - Pi) < angleErrer)
								{
									score += 1000;
								}
								//角度の合計が360度
								if (abs(ang2_1 + ang2_2 - Pi * 2) < angleErrer)
								{
									score += 1200;
								}
							}
						}
						//角度の合計が360度
						if (abs(ang1_1 + ang1_2 - Pi * 2) < angleErrer)
						{
							if (len2_1 - len2_2 + lengthErrer > 0)
							{
								if (ang2_1 < Pi + angleErrer)
								{
									score += 1200;
								}
							}
							else if (len2_2 - len2_1 + lengthErrer > 0)
							{
								if (ang2_2 < Pi + angleErrer)
								{
									score += 1200;
								}
							}
							else
							{
								//角度の合計が180度
								if (abs(ang2_1 + ang2_2 - Pi) < angleErrer)
								{
									score += 1600;
								}
								//角度の合計が360度
								if (abs(ang2_1 + ang2_2 - Pi * 2) < angleErrer)
								{
									score += 2000;
								}
							}
						}
					}
				}

				maxScore = max(maxScore, score);
			}
		}

		return maxScore;
	}

};

void Main()
{
	/*
	Image image(L"iphone.jpg");

	Window::Resize(image.size);

	image.grayscale().threshold(60);

	image.savePNG(L"二値化.png");

	image = Filter::filter(image);

	image.savePNG(L"ノイズ除去.png");
	*/

	Image image(L"ノイズ除去.png");

	image = Filter::median(image);

	auto images = ImageChip::chip(image);

	for (const auto& i : step(images.size())) images[i].savePNG(L"picture/分離" + Format(i) + L".png");

	Array<Array<Point>> pieces;

	ExtraPiece extra;

	for (auto& img : images)
	{
		const auto vec = extra.getPiece(img);
		for (size_t i = 0, size = vec.size(); i < size; i++)
			Line(vec[i], vec[(i + 1) % size]).overwrite(img, 4, Palette::Red);

		pieces.emplace_back(vec);
	}

	for (const auto& i : step(images.size()))
	{
		Window::Resize(images[i].size);
		System::Update();
		Texture tex(images[i]);
		tex.draw();
		ScreenCapture::Save(L"picture/形状取得" + Format(i) + L".png");
	}
	System::Update();

	AI ai;

	ai.think(pieces);

}
