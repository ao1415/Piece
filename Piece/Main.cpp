#include <Siv3D.hpp>

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

void Main()
{
	Image image(L"piece2.png");
	Texture texture(image);

	Window::Resize(texture.size * 2);

	ExtraPiece extra;

	while (System::Update())
	{
		texture.draw();

		const auto vec = extra.getPiece(image);

		for (size_t i = 0, size = vec.size(); i < size; i++)
		{
			const Point p1 = Point(texture.size.x, 0) + vec[i];
			const Point p2 = Point(texture.size.x, 0) + vec[(i + 1) % size];
			Line(p1, p2).draw(Palette::Yellow);
		}

		texture.draw(Point(texture.size.x / 2, texture.size.y), Alpha(127));
		for (size_t i = 0, size = vec.size(); i < size; i++)
		{
			const Point p1 = Point(texture.size.x / 2, texture.size.y) + vec[i];
			const Point p2 = Point(texture.size.x / 2, texture.size.y) + vec[(i + 1) % size];
			Line(p1, p2).draw(Palette::Yellow);
		}

	}

}
