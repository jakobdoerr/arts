#include "interactive_workspace.h"
#include "global_data.h"
#include "auto_workspace.h"
using global_data::md_data;

extern Verbosity verbosity_at_launch;
extern void (*getaways[])(Workspace&, const MRecord&);
extern WorkspaceMemoryHandler wsmh;
extern std::string string_buffer;

Index get_wsv_id(const char*);

size_t InteractiveWorkspace::n_anonymous_variables_ = 0;

InteractiveWorkspace::InteractiveWorkspace() : Workspace(), verbosity_(1, 1, 1)
{
    Workspace::initialize();
    verbosity_at_launch = verbosity_;
}

void InteractiveWorkspace::initialize() {
    define_wsv_group_names();
    Workspace::define_wsv_data();
    Workspace::define_wsv_map();
    define_md_data_raw();
    expand_md_data_raw_to_md_data();
    define_md_map();
    define_md_raw_map();
    define_agenda_data();
    define_agenda_map();
    assert( check_agenda_data() );
    define_species_data();
    define_species_map();
    define_lineshape_data();
    define_lineshape_norm_data();
}

const char * InteractiveWorkspace::execute_agenda(const Agenda *a)
{
    resize();
    try {
        a->execute(*this);
    } catch(const std::runtime_error &e) {
        string_buffer = e.what();
        return string_buffer.c_str();
    }
    return nullptr;
}

const char * InteractiveWorkspace::execute_workspace_method(long id,
                                                            const ArrayOfIndex &output,
                                                            const ArrayOfIndex &input)
{
    // Make sure verbosity is set.
    Index wsv_id_verbosity = get_wsv_id("verbosity");
    Verbosity& verbosity = *((Verbosity*)this->operator[](wsv_id_verbosity));
    verbosity.set_main_agenda(true);

    CREATE_OUTS;

    const MdRecord &m = md_data[id];
    if (m.SetMethod()) {
        swap(output[0], input[0]);
        return nullptr;
    }

    TokVal t{};
    Agenda a{};
    try {
        MRecord mr(id, output, input, t, a);
        getaways[id](*this, mr);
    } catch (const std::runtime_error &e) {
        string_buffer = e.what();
        return string_buffer.c_str();
    }
    return nullptr;
}

void InteractiveWorkspace::set_agenda_variable(Index id, const Agenda &src)
{
    *reinterpret_cast<Agenda*>(this->operator[](id)) = src;
    resize();
}

void InteractiveWorkspace::set_index_variable(Index id, const Index &src)
{
    *reinterpret_cast<Index*>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_numeric_variable(Index id, const Numeric &src)
{
    *reinterpret_cast<Numeric*>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_string_variable(Index id, const char *src)
{
    *reinterpret_cast<String*>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_array_of_string_variable(Index id,
                                                        size_t n,
                                                        const char * const *src)
{
    ArrayOfString *dst = reinterpret_cast<ArrayOfString*>(this->operator[](id));
    dst->resize(n);
    for (size_t i = 0; i < n; ++i) {
        dst->operator[](i) = String(src[i]);
    }
}

void InteractiveWorkspace::set_array_of_index_variable(Index id,
                                                       size_t n,
                                                       const Index *src)
{
    ArrayOfIndex *dst = reinterpret_cast<ArrayOfIndex*>(this->operator[](id));
    dst->resize(n);
    for (size_t i = 0; i < n; ++i) {
        dst->operator[](i) = src[i];
    }
}

void InteractiveWorkspace::set_vector_variable(Index id,
                                               size_t n,
                                               const Numeric *src)
{
    Vector *dst = reinterpret_cast<Vector*>(this->operator[](id));
    dst->resize(n);
    for (size_t i = 0; i < n; ++i) {
        dst->operator[](i) = src[i];
    }
}

void InteractiveWorkspace::set_matrix_variable(Index id,
                                               size_t m,
                                               size_t n,
                                               const Numeric *src)
{
    Matrix *dst = reinterpret_cast<Matrix*>(this->operator[](id));
    dst->resize(m, n);
    for (size_t i = 0; i < n * m; ++i) {
        dst->get_c_array()[i] = src[i];
    }
}

void InteractiveWorkspace::set_tensor3_variable(Index id,
                                                size_t l,
                                                size_t m,
                                                size_t n,
                                                const Numeric *src)
{
    Tensor3 *dst = reinterpret_cast<Tensor3*>(this->operator[](id));
    dst->resize(l, m, n);
    for (size_t i = 0; i < l * n * m; ++i) {
        dst->get_c_array()[i] = src[i];
    }
}

void InteractiveWorkspace::set_tensor4_variable(Index id,
                                                size_t k,
                                                size_t l,
                                                size_t m,
                                                size_t n,
                                                const Numeric *src)
{
  Tensor4 *dst = reinterpret_cast<Tensor4*>(this->operator[](id));
  dst->resize(k, l, m, n);
  for (size_t i = 0; i < k * l * n * m; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::resize()
{
    Array<stack<WsvStruct *>> ws_new(wsv_data.nelem());
    std::copy(ws.begin(), ws.end(), ws_new.begin());
    //std::copy(ws.end() - n_anonymous_variables_, ws.end(), ws_new.begin() + wsv_data.nelem());
    std::swap(ws, ws_new);
}

Index InteractiveWorkspace::add_variable(Index group_id, const char *name)
{
    if (wsv_data.size() != ws.size()) {
        resize();
    }
    Index id = static_cast<Index>(ws.size());

    ws.push_back(stack<WsvStruct *>());
    push(ws.size()-1, nullptr);
    ws.back().top()->wsv = wsmh.allocate (group_id);
    ws.back().top()->auto_allocated = true;
    ws.back().top()->initialized = true;

    String s;
    if (name) {
        s = String(name);
    } else {
        std::stringstream stream;
        stream << "anonymous_variable_" << n_anonymous_variables_;
        s = stream.str();
    }

    wsv_data.push_back(WsvRecord(s.c_str(), "Created by C API.", group_id));
    WsvMap[s] = id;

    ++n_anonymous_variables_;
    return id;
}

void InteractiveWorkspace::erase_variable(Index i, Index group_id)
{
    WsvStruct *wsvs;
    while (ws[i].size())
    {
        wsvs = ws[i].top();
        if (wsvs->auto_allocated && wsvs->wsv)
        {
            wsmh.deallocate (group_id, wsvs->wsv);
        }
        delete (wsvs);
        ws[i].pop ();
    }
    ws.erase(ws.begin() + i);

    WsvMap.erase(wsv_data[i].Name());
    wsv_data.erase(wsv_data.begin() + i);
    --n_anonymous_variables_;
}

void InteractiveWorkspace::swap(Index i, Index j)
{
    if (is_initialized(i) && is_initialized(j)) {
        std::swap(ws[i], ws[j]);
    } else if (!is_initialized(i) && is_initialized(j)) {
        ws[i].push(ws[j].top());
        ws[j].pop();
    } else if (is_initialized(i) && !is_initialized(j)) {
        ws[j].push(ws[i].top());
        ws[i].pop();
    }
}
