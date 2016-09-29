#include <Siv3D.hpp>
#include <stack>

using namespace std;

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


void Main()
{

	Webcam webcam;

	if (!webcam.open(1, Size(1920, 1080)))
		return;
	if (!webcam.start())
		return;

	Image original;
	while (System::Update())
	{
		if (webcam.hasNewFrame())
		{
			webcam.getFrame(original);
			break;
		}
	}

	if (!FileSystem::Exists(L"./chip"))
		FileSystem::CreateDirectories(L"./chip");

	original.savePNG(L"./chip/カメラ.png");

	Image extraction = original.grayscaled().thresholded(151);
	const auto images = ImageChip::chip(extraction);

	for (const auto& num : step(images.size()))
		images[num].savePNG(L"./chip/" + Pad(num, { L'0',2 }) + L".png");

	Console::Open();
	cout << images.size() << endl;
	system("pause");

}
