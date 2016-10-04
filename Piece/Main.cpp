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
			if (range > 15)
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
					const Point p = que.front();
					que.pop();

					const Point p1 = p + Point(0, -1);
					const Point p2 = p + Point(-1, 0);
					const Point p3 = p + Point(1, 0);
					const Point p4 = p + Point(0, 1);

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

				if (points.size() > 20)
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

};

void Main()
{
	Image image(L"抽出.png");

	ExtraPiece extra;

	image = Filter::filter(image);

	auto images = ImageChip::chip(image);

	for (const int i : step(images.size()))
	{
		System::Update();
		images[i].savePNG(L"picture/分離" + Format(i) + L".png");
	}

	return;

	for (auto& img : images)
	{
		const auto vec = extra.getPiece(img);

		for (size_t i = 0, size = vec.size(); i < size; i++)
		{
			Line(vec[i], vec[(i + 1) % size]).overwrite(img, Palette::Yellow);
		}
	}

	for (const int i : step(images.size()))
	{
		System::Update();
		Texture tex(images[i]);
		tex.draw();
		ScreenCapture::Save(L"picture/形状取得" + Format(i) + L".png");
	}

	while (System::Update())
	{

	}

}
