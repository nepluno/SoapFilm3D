// YImage.h
//

#ifdef HAVE_PNG

#ifndef __YImage_h__
#define __YImage_h__

// file loading/saving automatically picks up changes to this struct.
// The possibilities are: ARGB, ABGR, RGBA, BGRA.

struct YPixel {
  unsigned char r;
  unsigned char g;
  unsigned char b;
  unsigned char a;
};

class YImage {
 public:
  YImage();

  YImage(const YImage&);

  virtual ~YImage();

  YImage& operator=(const YImage&);

  bool save(const char* fname) const;

  bool load(const char* fname);

  YPixel* data();

  const YPixel* data() const;

  YPixel& at(int i, int j);

  const YPixel& at(int i, int j) const;

  int width() const;

  int height() const;

  void resize(int width, int height);

  // flip vertically
  void flip();

  // flip horizontally
  void mirror();

  // average rgb
  void greyscale();

 protected:
  int m_width;
  int m_height;
  YPixel* m_data;  // raw image data
};

#endif /* __YImage_h__ */

#endif
