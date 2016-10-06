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

	Array<std::pair<size_t, size_t>> think(const Array<Array<Point>>& points) {

		Array<std::pair<size_t, size_t>> pattern;

		getData(points);

		const size_t size = pieces.size();

		struct PairMatch {
			PairMatch() = default;
			PairMatch(size_t _n1, size_t _n2, int s) : n1(_n1), n2(_n2), score(s) {}

			size_t n1;
			size_t n2;
			int score;

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
			pattern.push_back({ data.n1,data.n2 });
			scoreQue.pop();
			printf("%3d, %3d: %d\n", (int)data.n1, (int)data.n2, data.score);
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

	const int match(const size_t n1, const size_t n2) const {

		const double angleErrer = 3.0 / 180 * Pi;
		const double lengthErrer = 10;

		const auto size1 = pieces[n1].corners.size();
		const auto size2 = pieces[n2].corners.size();

		int maxScore = 0;

		for (const auto i : step(size1))
		{
			for (const auto j : step(size2))
			{
				int score = eval(n1, n2, i, j);

				maxScore = max(maxScore, score);
			}
		}

		return maxScore;
	}

	const int eval(const size_t pieceId1, const size_t pieceId2, const size_t cornerId1, const size_t cornerId2) const {

		int score1 = 0;
		int score2 = 0;

		const double errerLenght = 10;
		const double errerAngle = Radians(3);

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

		flag = false;
		if (abs(len1_1 - len2_1) < errerLenght) score1 += 1;
		if (abs(len1_2 - len2_2) < errerLenght) { score1 += 3; flag = true; }
		if (abs(len1_3 - len2_3) < errerLenght) score1 += 1;

		if (abs(ang1_1 + ang2_1 - Pi) < errerAngle) score1 += 3;
		if (abs(ang1_2 + ang2_2 - Pi) < errerAngle) score1 += (flag ? 5 : 3);
		if (abs(ang1_3 + ang2_3 - Pi) < errerAngle) score1 += (flag ? 5 : 3);

		if (abs(ang1_1 + ang2_1 - Pi * 2) < errerAngle) score1 += 3;
		if (abs(ang1_2 + ang2_2 - Pi * 2) < errerAngle) score1 += (flag ? 5 : 3);
		if (abs(ang1_3 + ang2_3 - Pi * 2) < errerAngle) score1 += (flag ? 5 : 3);

		flag = false;
		if (abs(len1_1 - len2_3) < errerLenght) score2 += 1;
		if (abs(len1_2 - len2_2) < errerLenght) { score2 += 3; flag = true; }
		if (abs(len1_3 - len2_1) < errerLenght) score2 += 1;

		if (abs(ang1_1 + ang2_3 - Pi) < errerAngle) score2 += 3;
		if (abs(ang1_2 + ang2_2 - Pi) < errerAngle) score2 += (flag ? 5 : 3);
		if (abs(ang1_3 + ang2_1 - Pi) < errerAngle) score2 += (flag ? 5 : 3);

		if (abs(ang1_1 + ang2_3 - Pi * 2) < errerAngle) score2 += 3;
		if (abs(ang1_2 + ang2_2 - Pi * 2) < errerAngle) score2 += (flag ? 5 : 3);
		if (abs(ang1_3 + ang2_1 - Pi * 2) < errerAngle) score2 += (flag ? 5 : 3);

		return max(score1, score2);
	}

};

void Main()
{
	/*
	//Image image(L"iphone.jpg");
	//Window::Resize(image.size);

	Image image = Dialog::OpenImage();

	GUI gui(GUIStyle::Default);
	gui.add(L"sl", GUISlider::Create(0, 255, 60, true));
	gui.add(L"bu", GUIButton::Create(L"OK"));
	gui.show(true);

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

	const auto pattern = ai.think(pieces);
	//*/

	Window::Resize(1280, 960);
	while (System::Update())
	{
		const int W = 6, H = 4;
		for (int y = 0; y < H; y++)
		{
			for (int x = 0; x < W; x++)
			{
				Rect(Point(x * 1280 / W, y * 960 / H) + Point(10, 10), Point(1280 / W, 960 / H) - Point(20, 20)).drawFrame(1, 1, Palette::White);
			}
			Line({ 0,y * 960 / H }, { 320,y * 960 / H }).drawArrow(2);
			Line({ 320,y * 960 / H }, { 640,y * 960 / H }).drawArrow(2);
			Line({ 640,y * 960 / H }, { 960,y * 960 / H }).drawArrow(2);
			Line({ 960,y * 960 / H }, { 1280,y * 960 / H }).drawArrow(2);
		}
	}

}
