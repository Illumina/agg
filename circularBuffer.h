class circularBuffer {
 public:

  circularBuffer() {
    resize(0);
  }

  circularBuffer(int buffer_size) {
    resize(buffer_size);
  }

  void resize(int buffer_size) {
    _bufsize=buffer_size;
    _offset=_bufn=0;
    _buf = new depthInterval[_bufsize];
    //_buf.resize(_bufsize);
  }

  ~circularBuffer() {    delete[] _buf;  }

  int size() {    return(_bufn);  }
  bool  empty() {    return(_bufn==0);  }
  bool  full() {    return(_bufn==_bufsize);  }

  bool push_back(depthInterval o) {
    int idx=(_offset + _bufn)%_bufsize;
    //    cerr << "push_back " << _offset << " " <<_bufn<<" "<<idx<<" "<<endl;
    assert(idx<_bufsize&&idx>=0);
    if(_bufn<_bufsize) _buf[idx]=o;
    else die("buffer overflow");
    _bufn++;
    return(true);
  }

  depthInterval *front() { 
    assert(_bufn>0);
    assert(_offset>=0 && _offset<_bufsize);      
    return(&_buf[_offset]);
  }

  depthInterval *back() { 
    assert(_bufn>0);
    return(&_buf[(_offset+_bufn-1)%_bufsize]);
  }

  void pop_front() {
    if(_bufn>0) {
      _offset++;
      _offset %= _bufsize;
      _bufn--;
    }
    else 
      die("buffer empty (pop_front)");
  }
 private:
  //  vector<depthInterval> _buf;
  depthInterval *_buf;
  int _bufsize,_offset,_bufn;  
};


