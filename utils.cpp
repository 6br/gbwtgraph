#include <gbwtgraph/utils.h>

#include <sstream>
#include <vector>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Other class variables.

const std::string SequenceSource::TRANSLATION_EXTENSION = ".trans";

//------------------------------------------------------------------------------

std::string
Version::str(bool verbose)
{
  std::ostringstream ss;
  if(verbose) { ss << "GBWTGraph version "; }
  else { ss << "v"; }
  ss << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_VERSION;
  if(verbose) { ss << " (file format version " << GRAPH_VERSION << ")"; }
  return ss.str();
}

void
Version::print(std::ostream& out, const std::string& tool_name, bool verbose, size_t new_lines)
{
  out << tool_name;
  if(verbose) { out << std::endl; }
  else { out << " "; }
  out << str(verbose);
  for(size_t i = 0; i < new_lines; i++) { out << std::endl; }
}

//------------------------------------------------------------------------------

const std::vector<char> complement =
{
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', '$', '#', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', '-', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

  'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C',   'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N',
  'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S',   'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 't', 'v', 'g', 'h', 'N', 'N', 'c',   'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N',
  'N', 'q', 'y', 'w', 'a', 'a', 'b', 's',   'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N',

  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

std::string
reverse_complement(const std::string& seq)
{
  std::string result = seq;
  reverse_complement_in_place(result);
  return result;
}

void
reverse_complement_in_place(std::string& seq)
{
  size_t swap_size = seq.size() / 2;
  for(size_t i = 0, j = seq.size() - 1; i < swap_size; i++, j--)
  {
    char tmp = seq[i];
    seq[i] = complement[seq[j]];
    seq[j] = complement[tmp];
  }
  
  if(seq.size() % 2 != 0)
  {
    seq[swap_size] = complement[seq[swap_size]];
  }
}

//------------------------------------------------------------------------------

void
SequenceSource::swap(SequenceSource& another)
{
  if(&another == this) { return; }

  this->nodes.swap(another.nodes);
  this->sequences.swap(another.sequences);
  this->segment_translation.swap(another.segment_translation);
  std::swap(this->next_id, another.next_id);
}

void
SequenceSource::add_node(nid_t node_id, const std::string& sequence)
{
  if(sequence.empty()) { return; }
  handle_t handle = this->get_handle(node_id, false);
  if(this->nodes.find(handle) != this->nodes.end()) { return; }
  size_t offset = this->sequences.size();
  this->sequences.insert(this->sequences.end(), sequence.begin(), sequence.end());
  this->nodes[handle] = std::pair<size_t, size_t>(offset, sequence.length());
}

void
SequenceSource::add_node(nid_t node_id, view_type sequence)
{
  if(sequence.second == 0) { return; }
  handle_t handle = this->get_handle(node_id, false);
  if(this->nodes.find(handle) != this->nodes.end()) { return; }
  size_t offset = this->sequences.size();
  this->sequences.insert(this->sequences.end(), sequence.first, sequence.first + sequence.second);
  this->nodes[handle] = std::pair<size_t, size_t>(offset, sequence.second);
}

void
SequenceSource::translate_segment(const std::string& name, view_type sequence, size_t max_length)
{
  if(this->segment_translation.find(name) != this->segment_translation.end() || sequence.second == 0) { return; }

  nid_t start = this->next_id;
  nid_t limit = start + (sequence.second + max_length - 1) / max_length;
  for(nid_t id = start; id < limit; id++)
  {
    size_t offset = (id - start) * max_length;
    size_t length = std::min(max_length, sequence.second - offset);
    this->add_node(id, view_type(sequence.first + offset, length));
  }

  this->segment_translation[name] = std::make_pair(start, limit);
  this->next_id = limit;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph

