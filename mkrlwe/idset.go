package mkrlwe

type IDSet struct {
	Value map[string]struct{} //empty structs occupy 0 memory
}

func (s *IDSet) Has(v string) bool {
	_, ok := s.Value[v]
	return ok
}

func (s *IDSet) Add(v string) {
	if v == "0" {
		panic("Cannot IDSet Add : 0 cannot be used")
	}
	s.Value[v] = struct{}{}
}

func (s *IDSet) Remove(v string) {
	delete(s.Value, v)
}

func (s *IDSet) Size() int {
	return len(s.Value)
}

func NewIDSet() *IDSet {
	s := &IDSet{}
	s.Value = make(map[string]struct{})
	return s
}

func (s *IDSet) CopyNew() *IDSet {
	res := NewIDSet()
	for v := range s.Value {
		res.Add(v)
	}
	return res
}

func (s *IDSet) Union(s2 *IDSet) *IDSet {
	res := NewIDSet()
	for v := range s.Value {
		res.Add(v)
	}

	for v := range s2.Value {
		res.Add(v)
	}
	return res
}
