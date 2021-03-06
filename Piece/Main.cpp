﻿#include <Siv3D.hpp>
#include <stack>
#include <queue>

using namespace std;

class ExtraPiece {
public:

	const Array<Point> getPiece(const Image& image, const double max_distance = 3.0) {

		const double errerAngle = Radians(15);

		const Polygon polygon = Imaging::FindContour(image, false);
		const Array<Vec2> vecArray = polygon.simplified(max_distance).outer();

		size_t vecSize;
		Array<Vec2> vec1 = vecArray;
		Array<Vec2> vec2;

		std::cout << "元データ\t" << vec1.size() << std::endl;
		for (const auto& v : vec1) std::cout << v << std::endl;

		do
		{
			vecSize = vec1.size();
			for (size_t i = 0, size = vec1.size(); i < size; i++)
			{
				const Vec2 p1 = vec1[(i + 0) % size];
				const Vec2 p2 = vec1[(i + 1) % size];
				const Vec2 p3 = vec1[(i + 2) % size];

				const double ang1 = Math::Atan2(p2.y - p1.y, p2.x - p1.x);
				const double ang2 = Math::Atan2(p3.y - p2.y, p3.x - p2.x);

				if (abs(ang1 - ang2) > errerAngle || abs(ang1 + ang2 - Pi) > errerAngle)
					vec2.push_back(Point(p1.x, p1.y));
				else
				{
					vec2.push_back(Point(p1.x, p1.y));
					i++;
				}
			}
			vec1.swap(vec2);
			vec2.clear();

		} while (vecSize > vec1.size());

		std::cout << "角度調節1\t" << vec1.size() << std::endl;
		for (const auto& v : vec1) std::cout << v << std::endl;

		do
		{
			vecSize = vec1.size();
			for (size_t i = 0, size = vec1.size(); i < size; i++)
			{
				const double range = Math::Sqrt(Math::Square(vec1[i].x - vec1[(i + 1) % size].x) + Math::Square(vec1[i].y - vec1[(i + 1) % size].y));
				const double angel = Math::Atan2(vec1[(i + 1) % size].y - vec1[i].y, vec1[(i + 1) % size].x - vec1[i].x);


				//最低の辺の長さ未満ならば辺を削除する
				if (range > 15)
				{
					vec2.push_back(vec1[i]);
				}
				else
				{
					const Vec2 p1 = vec1[(i + size - 1) % size];
					const Vec2 p2 = vec1[(i + 0) % size];
					const Vec2 p3 = vec1[(i + 1) % size];
					const Vec2 p4 = vec1[(i + 2) % size];

					const Vec2 v1 = p2 - p1;
					const Vec2 v2 = p3 - p4;

					const double sinv1 = (v1.x != 0) ? (v1.y / v1.x) : (1);
					const double sinv2 = (v2.x != 0) ? (v2.y / v2.x) : (1);

					const double x1 = ((p1.x*sinv1 - p4.x*sinv2) - (p1.y - p4.y)) / (sinv1 - sinv2);
					const double y1 = p1.y + sinv1*(x1 - p1.x);
					const double x2 = ((p1.x*sinv1 - p4.x*sinv2) - (p1.y - p4.y)) / (sinv1 - sinv2);
					const double y2 = p4.y + sinv2*(x2 - p4.x);

					i++;

					vec2.emplace_back((int)x1, (int)y1);
				}
			}

			vec1.swap(vec2);
			vec2.clear();

		} while (vecSize > vec1.size());

		std::cout << "辺調節\t" << vec1.size() << std::endl;
		for (const auto& v : vec1) std::cout << v << std::endl;


		do
		{
			vecSize = vec1.size();
			for (size_t i = 0, size = vec1.size(); i < size; i++)
			{
				const Vec2 p1 = vec1[(i + 0) % size];
				const Vec2 p2 = vec1[(i + 1) % size];
				const Vec2 p3 = vec1[(i + 2) % size];

				const double ang1 = Math::Atan2(p2.y - p1.y, p2.x - p1.x);
				const double ang2 = Math::Atan2(p3.y - p2.y, p3.x - p2.x);

				if (abs(ang1 - ang2) > errerAngle || abs(ang1 + ang2 - Pi) > errerAngle)
					vec2.push_back(Point(p1.x, p1.y));
				else
				{
					vec2.push_back(Point(p1.x, p1.y));
					i++;
				}
			}
			vec1.swap(vec2);
			vec2.clear();

		} while (vecSize > vec1.size());

		std::cout << "角度調節2\t" << vec1.size() << std::endl;
		for (const auto& i : step(vec1.size())) std::cout << vec1[i] << "\t," << Math::Atan2(vec1[(i + 1) % vec1.size()].y - vec1[(i + 0) % vec1.size()].y, vec1[(i + 1) % vec1.size()].x - vec1[(i + 0) % vec1.size()].x) << std::endl;


		Array<Point> points;
		for (const auto& v : vec1)
			points.emplace_back(Point(v.x, v.y));

		return points;
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

				if (points.size() < 1000)
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

struct PairMatch {
	PairMatch() = default;
	PairMatch(size_t _n1, size_t _n2, size_t _v1, size_t _v2, int s, double a, Point p, bool r) : n1(_n1), n2(_n2), v1(_v1), v2(_v2), score(s), angle(a), pos(p), reverce(r) {}

	size_t n1;
	size_t n2;
	size_t v1;
	size_t v2;
	int score;

	//もう一つのピース
	double angle;
	Point pos;
	bool reverce;

	const bool operator<(const PairMatch& other) const {
		if (score < other.score) return true;
		else if (score > other.score) return false;
		if (n1 < other.n1) return false;
		else if (n1 > other.n1) return true;
		if (n2 < other.n2) return false;
		else if (n2 > other.n2) return true;
		return false;
	}

};

class AI {
public:

	Array<PairMatch> think(const Array<Array<Point>>& points) {

		Array<PairMatch> pattern;

		getData(points);

		const size_t size = pieces.size();

		priority_queue<PairMatch> scoreQue;
		vector<PairMatch> scoreVec;

		for (size_t i = 0; i < size; i++)
		{
			for (size_t j = i + 1; j < size; j++)
			{
				pair<int, int> pair1;
				pair<int, int> pair2;

				const int score1 = match(i, j, pair1);
				const int score2 = 0;//match(i, j, true, pair2);
				const int score = max(score1, score2);

				if (score > 0)
				{
					if (score1 > score2)
					{
						const size_t size1 = pieces[i].corners.size();
						const size_t size2 = pieces[j].corners.size();

						const Point pos1_1 = pieces[i].corners[(pair1.first + 0) % size1].point;
						const Point pos1_2 = pieces[i].corners[(pair1.first + 1) % size1].point;
						const Point pos2_1 = pieces[j].corners[(pair1.second + 0) % size2].point;
						const Point pos2_2 = pieces[j].corners[(pair1.second + 1) % size2].point;

						const double ang1 = Math::Atan2((double)pos1_2.y - pos1_1.y, (double)pos1_1.x - pos1_1.x);
						const double ang2 = Math::Atan2((double)pos2_2.y - pos2_1.y, (double)pos2_1.x - pos2_1.x);

						const double dif = ang1 - ang2;

						const Point rot = Point(pos2_1.x*cos(dif) - pos2_1.y*sin(dif), pos2_1.x*sin(dif) + pos2_1.y*cos(dif));
						const Point difP = pos1_1 - rot;

						scoreVec.push_back(PairMatch(i, j, pair1.first, pair1.second, score, dif, difP, false));
						scoreQue.push(PairMatch(i, j, pair1.first, pair1.second, score, dif, difP, false));
					}
					else
					{
						const size_t size1 = pieces[i].corners.size();
						const size_t size2 = pieces[j].corners.size();

						const Point pos1_1 = pieces[i].corners[(pair1.first + 0) % size1].point;
						const Point pos1_2 = pieces[i].corners[(pair1.first + 1) % size1].point;
						const Point pos2_1 = pieces[j].corners[(pair1.second + 2) % size2].point;
						const Point pos2_2 = pieces[j].corners[(pair1.second + 1) % size2].point;

						const double ang1 = Math::Atan2((double)pos1_2.y - pos1_1.y, (double)pos1_1.x - pos1_1.x);
						const double ang2 = Math::Atan2((double)pos2_2.y - pos2_1.y, (double)pos2_1.x - pos2_1.x);

						const double dif = ang1 - ang2;

						const Point rot = Point(pos2_1.x*cos(dif) - pos2_1.y*sin(dif), pos2_1.x*sin(dif) + pos2_1.y*cos(dif));
						const Point difP = pos1_1 - rot;

						scoreVec.push_back(PairMatch(i, j, pair2.first, pair2.second, score, dif, difP, true));
						scoreQue.push(PairMatch(i, j, pair2.first, pair2.second, score, dif, difP, true));
					}

				}


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
			pattern.emplace_back(data);
		}

		return pattern;
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

	const int match(const size_t n1, const size_t n2, pair<int, int>& maxRad) const {

		const auto size1 = pieces[n1].corners.size();
		const auto size2 = pieces[n2].corners.size();

		int maxScore = 0;

		for (const auto i : step(size1))
		{
			for (const auto j : step(size2))
			{
				int score = eval(n1, n2, i, j);

				if (maxScore < score)
				{
					maxScore = score;
					maxRad = { (int)i,(int)j };
				}
			}
		}

		return maxScore;
	}

	const int eval(const size_t pieceId1, const size_t pieceId2, const size_t cornerId1, const size_t cornerId2) const {

		int score1 = 0;
		int score2 = 0;

		const double errerLenght = 10;
		const double errerAngle = Radians(10);

		const auto size1 = pieces[pieceId1].corners.size();
		const auto size2 = pieces[pieceId2].corners.size();

		const double len1_1 = pieces[pieceId1].corners[(cornerId1 + 0) % size1].length;
		const double len1_2 = pieces[pieceId1].corners[(cornerId1 + 1) % size1].length;
		const double len1_3 = pieces[pieceId1].corners[(cornerId1 + 2) % size1].length;

		const double len2_1 = pieces[pieceId2].corners[(cornerId2 + 0) % size2].length;
		const double len2_2 = pieces[pieceId2].corners[(cornerId2 + 1) % size2].length;
		const double len2_3 = pieces[pieceId2].corners[(cornerId2 + 2) % size2].length;

		const double ang1_1 = pieces[pieceId1].corners[(cornerId1 + 0) % size1].angle;
		const double ang1_2 = pieces[pieceId1].corners[(cornerId1 + 1) % size1].angle;
		const double ang1_3 = pieces[pieceId1].corners[(cornerId1 + 2) % size1].angle;

		const double ang2_1 = pieces[pieceId2].corners[(cornerId2 + 0) % size2].angle;
		const double ang2_2 = pieces[pieceId2].corners[(cornerId2 + 1) % size2].angle;
		const double ang2_3 = pieces[pieceId2].corners[(cornerId2 + 2) % size2].angle;

		bool flag;

		if (!false)
		{
			flag = false;
			/*
			if (abs(len1_1 - len2_3) < errerLenght) score1 += 100;
			if (abs(len1_2 - len2_2) < errerLenght) { score1 += 150; flag = true; }
			if (abs(len1_3 - len2_1) < errerLenght) score1 += 100;

			if (abs(ang1_1 + ang2_3 - Pi) < errerAngle) score1 += 100;
			if (abs(ang1_2 + ang2_2 - Pi) < errerAngle) score1 += (flag ? 150 : 100);
			if (abs(ang1_3 + ang2_1 - Pi) < errerAngle) score1 += (flag ? 150 : 100);

			if (abs(ang1_1 + ang2_3 - Pi * 2) < errerAngle) score1 += 100;
			if (abs(ang1_2 + ang2_2 - Pi * 2) < errerAngle) score1 += (flag ? 300 : 150);
			if (abs(ang1_3 + ang2_1 - Pi * 2) < errerAngle) score1 += (flag ? 300 : 150);
			*/

			if (abs(len1_1 - len2_3) < errerLenght)
				if (abs(len1_2 - len2_2) < errerLenght)
					if (abs(len1_3 - len2_1) < errerLenght)
						score1 += 50;

			if (abs(ang1_1 + ang2_3 - Pi * 2) < errerAngle)
			{
				if (abs(ang1_2 + ang2_2 - Pi * 2) < errerAngle)
				{
					if (abs(ang1_3 + ang2_1 - Pi * 2) < errerAngle)
						score1 += 100;
					else
						score1 += 50;
				}
			}
			else if (abs(ang1_2 + ang2_2 - Pi * 2) < errerAngle)
			{
				if (abs(ang1_3 + ang2_1 - Pi * 2) < errerAngle)
					score1 += 50;
			}

			if (abs(ang1_2 + ang2_2 - Pi * 2) < errerAngle)
			{
				if (abs(len1_2 - len2_2) < errerLenght)
				{
					if (abs(len1_1 - len2_3) < errerLenght)
						score1 += 75;
					else
						score1 += 25;
				}
			}
			if (abs(ang1_3 + ang2_1 - Pi * 2) < errerAngle)
			{
				if (abs(len1_2 - len2_2) < errerLenght)
				{
					if (abs(len1_3 - len2_1) < errerLenght)
						score1 += 75;
					else
						score1 += 25;
				}
			}

		}

		return max(score1, score2);
	}

};

void Main()
{
	//*
	//Image image(L"iphone.jpg");
	//Window::Resize(image.size);

	Image image = Dialog::OpenImage();

	GUI gui(GUIStyle::Default);
	gui.add(L"sl", GUISlider::Create(0, 255, 60, true));
	gui.add(L"bu", GUIButton::Create(L"OK"));
	gui.show(true);

	Window::Resize(1920, 1080);
	Graphics::SetBackground(Palette::Lightgreen);

	Texture texture = Texture(image.grayscaled().thresholded(60));
	while (System::Update())
	{
		if (gui.slider(L"sl").hasChanged)
		{
			texture = Texture(image.grayscaled().thresholded((uint8)gui.slider(L"sl").value));
		}

		if (gui.button(L"bu").pushed)
			break;

		texture.resize(Window::Size()).draw();
	}
	gui.show(false);

	image.grayscale().threshold((uint8)gui.slider(L"sl").value);

	image.savePNG(L"二値化.png");

	image = Filter::filter(image);

	image.savePNG(L"ノイズ除去.png");

	//Image image(L"ノイズ除去.png");

	image = Filter::median(image);

	auto images = ImageChip::chip(image);

	for (const auto& i : step(images.size()))
	{
		images[i].savePNG(L"picture/分離" + Format(i) + L".png");
	}

	Array<Array<Point>> pieces;

	ExtraPiece extra;

	for (const auto& c : step(images.size()))
	{
		std::cout << "ピース:" << c << std::endl;
		const auto vec = extra.getPiece(images[c]);
		for (size_t i = 0, size = vec.size(); i < size; i++)
			Line(vec[i], vec[(i + 1) % size]).overwriteArrow(images[c], 3, { 10.0,10.0 }, Palette::Red);

		pieces.emplace_back(vec);
	}

	/*
	for (const auto& i : step(images.size()))
	{
		Window::Resize(images[i].size);
		System::Update();
		Texture tex(images[i]);
		tex.draw();
		ScreenCapture::Save(L"picture/形状取得" + Format(i) + L".png");
	}
	System::Update();
	*/

	AI ai;

	const auto pattern = ai.think(pieces);
	//*/

	int scroll = 0;
	Array<Texture> textures;

	for (const auto& i : step(images.size()))
	{
		images[i].forEach([](Color& c) {if (c.r < 127) c.a = 0; });
		textures.push_back(Texture(images[i]));
	}

	Font font(24);

	while (System::Update())
	{
		const int W = 3, H = 3;
		for (int y = 0; y < H; y++)
		{
			for (int x = 0; x < W; x++)
			{
				const auto size = Point(Window::Size().x / W, Window::Size().y / H) - Point(20, 20);
				const Size hSize = { size.x / 2,size.y };
				const auto pos = Point(x * Window::Size().x / W, y * Window::Size().y / H) + Point(10, 10);
				const Point center = size / 2;

				if ((y + scroll)*W + x < pattern.size())
				{
					const auto pieceId1 = pattern[(y + scroll)*W + x].n1;
					const auto pieceId2 = pattern[(y + scroll)*W + x].n2;

					const auto p1_1 = pieces[pieceId1][(pattern[(y + scroll)*W + x].v1 + 0) % pieces[pieceId1].size()];
					const auto p1_2 = pieces[pieceId1][(pattern[(y + scroll)*W + x].v1 + 1) % pieces[pieceId1].size()];
					const auto p1_3 = pieces[pieceId1][(pattern[(y + scroll)*W + x].v1 + 2) % pieces[pieceId1].size()];
					const auto p1_4 = pieces[pieceId1][(pattern[(y + scroll)*W + x].v1 + 3) % pieces[pieceId1].size()];

					const auto p2_1 = pieces[pieceId2][(pattern[(y + scroll)*W + x].v2 + 0) % pieces[pieceId2].size()];
					const auto p2_2 = pieces[pieceId2][(pattern[(y + scroll)*W + x].v2 + 1) % pieces[pieceId2].size()];
					const auto p2_3 = pieces[pieceId2][(pattern[(y + scroll)*W + x].v2 + 2) % pieces[pieceId2].size()];
					const auto p2_4 = pieces[pieceId2][(pattern[(y + scroll)*W + x].v2 + 3) % pieces[pieceId2].size()];

					const double scale = min(
						min(
							min((double)hSize.x / textures[pieceId1].size.x,
							(double)(hSize.y - 68) / textures[pieceId1].size.y),
							min((double)hSize.x / textures[pieceId2].size.x,
							(double)(hSize.y - 68) / textures[pieceId2].size.y)
						),
						1.0
					);

					textures[pieceId1].scale(scale).draw(pos + Point(0, 0) + Point(10, 10));
					textures[pieceId2].scale(scale).draw(pos + Point(center.x, 0) + Point(10, 10));

					Line(p1_1*scale + pos + Point(0, 0) + Point(10, 10), p1_2*scale + pos + Point(0, 0) + Point(10, 10)).drawArrow(4, { 10,10 }, Palette::Black);
					Line(p1_2*scale + pos + Point(0, 0) + Point(10, 10), p1_3*scale + pos + Point(0, 0) + Point(10, 10)).drawArrow(4, { 10,10 }, Palette::Black);
					Line(p1_3*scale + pos + Point(0, 0) + Point(10, 10), p1_4*scale + pos + Point(0, 0) + Point(10, 10)).drawArrow(4, { 10,10 }, Palette::Black);

					Line(p2_1*scale + pos + Point(center.x, 0) + Point(10, 10), p2_2*scale + pos + Point(center.x, 0) + Point(10, 10)).drawArrow(4, { 10,10 }, Palette::Black);
					Line(p2_2*scale + pos + Point(center.x, 0) + Point(10, 10), p2_3*scale + pos + Point(center.x, 0) + Point(10, 10)).drawArrow(4, { 10,10 }, Palette::Black);
					Line(p2_3*scale + pos + Point(center.x, 0) + Point(10, 10), p2_4*scale + pos + Point(center.x, 0) + Point(10, 10)).drawArrow(4, { 10,10 }, Palette::Black);

					//textures[pieceId2].rotate(pattern[(y + scroll)*W + x].angle).draw(pos + center - pattern[(y + scroll)*W + x].pos);

					font(Format(L"(", Pad(pieceId1, { L'0',2 }), L", ", Pad(pieceId2, { L'0',2 }), L"):", pattern[(y + scroll)*W + x].score, L"\t", (y + scroll)*W + x, L"/", pattern.size())).draw({ pos.x,pos.y + size.y - 48 }, Palette::Black);

					//Window::SetTitle(Format(Pad(pieceId1, { L'0',2 }), L"-", Pad(pieceId2, { L'0',2 }), L":", pattern[(y + scroll)*W + x].score, L"-", scroll, L"/", pattern.size() / W / H));

					Rect(pos, size).drawFrame(1, 1, Palette::White);
				}
			}
		}

		if (Input::KeyUp.clicked)
			scroll = max(scroll - H, 0);
		else if (Input::KeyDown.clicked)
			scroll = min(scroll + H, (int)pattern.size() / W);

	}

}
