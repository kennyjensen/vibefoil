#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xbl.h"
#include "xoper.h"

namespace {

struct Json {
    enum class Type { kNull, kBool, kNumber, kString, kArray, kObject };
    Type type = Type::kNull;
    bool boolean = false;
    double number = 0.0;
    std::string str;
    std::vector<Json> array;
    std::unordered_map<std::string, Json> object;
};

void skip_ws(const std::string &s, size_t &i) {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) {
        ++i;
    }
}

std::string parse_string(const std::string &s, size_t &i) {
    if (s[i] != '"') {
        throw std::runtime_error("Expected string");
    }
    ++i;
    std::string out;
    while (i < s.size()) {
        const char c = s[i++];
        if (c == '"') {
            return out;
        }
        if (c == '\\') {
            if (i >= s.size()) {
                throw std::runtime_error("Bad escape");
            }
            const char e = s[i++];
            switch (e) {
                case '"':
                    out.push_back('"');
                    break;
                case '\\':
                    out.push_back('\\');
                    break;
                case '/':
                    out.push_back('/');
                    break;
                case 'b':
                    out.push_back('\b');
                    break;
                case 'f':
                    out.push_back('\f');
                    break;
                case 'n':
                    out.push_back('\n');
                    break;
                case 'r':
                    out.push_back('\r');
                    break;
                case 't':
                    out.push_back('\t');
                    break;
                default:
                    throw std::runtime_error("Unsupported escape");
            }
            continue;
        }
        out.push_back(c);
    }
    throw std::runtime_error("Unterminated string");
}

Json parse_value(const std::string &s, size_t &i);

Json parse_array(const std::string &s, size_t &i) {
    Json out;
    out.type = Json::Type::kArray;
    ++i;
    skip_ws(s, i);
    if (i < s.size() && s[i] == ']') {
        ++i;
        return out;
    }
    while (i < s.size()) {
        out.array.push_back(parse_value(s, i));
        skip_ws(s, i);
        if (i >= s.size()) {
            throw std::runtime_error("Unterminated array");
        }
        if (s[i] == ',') {
            ++i;
            skip_ws(s, i);
            continue;
        }
        if (s[i] == ']') {
            ++i;
            return out;
        }
        throw std::runtime_error("Bad array separator");
    }
    throw std::runtime_error("Unterminated array");
}

Json parse_object(const std::string &s, size_t &i) {
    Json out;
    out.type = Json::Type::kObject;
    ++i;
    skip_ws(s, i);
    if (i < s.size() && s[i] == '}') {
        ++i;
        return out;
    }
    while (i < s.size()) {
        skip_ws(s, i);
        const std::string key = parse_string(s, i);
        skip_ws(s, i);
        if (i >= s.size() || s[i] != ':') {
            throw std::runtime_error("Expected ':'");
        }
        ++i;
        skip_ws(s, i);
        out.object.emplace(key, parse_value(s, i));
        skip_ws(s, i);
        if (i >= s.size()) {
            throw std::runtime_error("Unterminated object");
        }
        if (s[i] == ',') {
            ++i;
            skip_ws(s, i);
            continue;
        }
        if (s[i] == '}') {
            ++i;
            return out;
        }
        throw std::runtime_error("Bad object separator");
    }
    throw std::runtime_error("Unterminated object");
}

Json parse_number(const std::string &s, size_t &i) {
    const char *start = s.c_str() + i;
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(start, &end);
    if (start == end || errno != 0) {
        throw std::runtime_error("Bad number");
    }
    i += static_cast<size_t>(end - start);
    Json out;
    out.type = Json::Type::kNumber;
    out.number = value;
    return out;
}

Json parse_value(const std::string &s, size_t &i) {
    skip_ws(s, i);
    if (i >= s.size()) {
        throw std::runtime_error("Unexpected end of input");
    }
    const char c = s[i];
    if (c == '{') {
        return parse_object(s, i);
    }
    if (c == '[') {
        return parse_array(s, i);
    }
    if (c == '"') {
        Json out;
        out.type = Json::Type::kString;
        out.str = parse_string(s, i);
        return out;
    }
    if (c == 't' && s.compare(i, 4, "true") == 0) {
        i += 4;
        Json out;
        out.type = Json::Type::kBool;
        out.boolean = true;
        return out;
    }
    if (c == 'f' && s.compare(i, 5, "false") == 0) {
        i += 5;
        Json out;
        out.type = Json::Type::kBool;
        out.boolean = false;
        return out;
    }
    if (c == 'n' && s.compare(i, 4, "null") == 0) {
        i += 4;
        Json out;
        out.type = Json::Type::kNull;
        return out;
    }
    return parse_number(s, i);
}

Json parse_json(const std::string &s) {
    size_t i = 0;
    Json out = parse_value(s, i);
    skip_ws(s, i);
    if (i != s.size()) {
        throw std::runtime_error("Trailing data");
    }
    return out;
}

void require_type(const Json &v, Json::Type type, const std::string &name) {
    if (v.type != type) {
        throw std::runtime_error("Type mismatch for " + name);
    }
}

void assign_scalar(double &out, const Json &v) {
    if (v.type == Json::Type::kNumber) {
        out = v.number;
        return;
    }
    throw std::runtime_error("Expected number");
}

void assign_scalar(int &out, const Json &v) {
    if (v.type == Json::Type::kNumber) {
        out = static_cast<int>(v.number);
        return;
    }
    throw std::runtime_error("Expected number");
}

void assign_scalar(bool &out, const Json &v) {
    if (v.type == Json::Type::kBool) {
        out = v.boolean;
        return;
    }
    if (v.type == Json::Type::kNumber) {
        out = (v.number != 0.0);
        return;
    }
    throw std::runtime_error("Expected bool");
}

void assign_scalar(std::string &out, const Json &v) {
    if (v.type == Json::Type::kString) {
        out = v.str;
        return;
    }
    throw std::runtime_error("Expected string");
}

void assign_vector(std::vector<double> &out, const Json &v) {
    require_type(v, Json::Type::kArray, "vector<double>");
    if (out.size() != v.array.size()) {
        throw std::runtime_error("vector<double> length mismatch");
    }
    for (size_t i = 0; i < v.array.size(); ++i) {
        out[i] = v.array[i].number;
    }
}

void assign_vector(std::vector<int> &out, const Json &v) {
    require_type(v, Json::Type::kArray, "vector<int>");
    if (out.size() != v.array.size()) {
        throw std::runtime_error("vector<int> length mismatch");
    }
    for (size_t i = 0; i < v.array.size(); ++i) {
        out[i] = static_cast<int>(v.array[i].number);
    }
}

void assign_vector(std::vector<bool> &out, const Json &v) {
    require_type(v, Json::Type::kArray, "vector<bool>");
    if (out.size() != v.array.size()) {
        throw std::runtime_error("vector<bool> length mismatch");
    }
    for (size_t i = 0; i < v.array.size(); ++i) {
        const Json &el = v.array[i];
        out[i] = (el.type == Json::Type::kBool) ? el.boolean : (el.number != 0.0);
    }
}

void assign_vector(std::vector<std::vector<double>> &out, const Json &v) {
    require_type(v, Json::Type::kArray, "vector<vector<double>>");
    if (out.size() != v.array.size()) {
        throw std::runtime_error("vector<vector<double>> length mismatch");
    }
    for (size_t i = 0; i < v.array.size(); ++i) {
        assign_vector(out[i], v.array[i]);
    }
}

void assign_vector(std::vector<std::vector<int>> &out, const Json &v) {
    require_type(v, Json::Type::kArray, "vector<vector<int>>");
    if (out.size() != v.array.size()) {
        throw std::runtime_error("vector<vector<int>> length mismatch");
    }
    for (size_t i = 0; i < v.array.size(); ++i) {
        assign_vector(out[i], v.array[i]);
    }
}

void assign_vector(std::vector<std::vector<std::vector<double>>> &out, const Json &v) {
    require_type(v, Json::Type::kArray, "vector<vector<vector<double>>>");
    if (out.size() != v.array.size()) {
        throw std::runtime_error("vector<vector<vector<double>>> length mismatch");
    }
    for (size_t i = 0; i < v.array.size(); ++i) {
        assign_vector(out[i], v.array[i]);
    }
}

void write_json_string(std::ostream &out, const std::string &value) {
    out << '"';
    for (char c : value) {
        switch (c) {
            case '"':
                out << "\\\"";
                break;
            case '\\':
                out << "\\\\";
                break;
            case '\b':
                out << "\\b";
                break;
            case '\f':
                out << "\\f";
                break;
            case '\n':
                out << "\\n";
                break;
            case '\r':
                out << "\\r";
                break;
            case '\t':
                out << "\\t";
                break;
            default:
                out << c;
        }
    }
    out << '"';
}

void dump_vector(std::ostream &out, const std::vector<double> &value) {
    out << '[';
    for (size_t i = 0; i < value.size(); ++i) {
        if (i) out << ',';
        out << std::setprecision(17) << value[i];
    }
    out << ']';
}

void dump_vector(std::ostream &out, const std::vector<int> &value) {
    out << '[';
    for (size_t i = 0; i < value.size(); ++i) {
        if (i) out << ',';
        out << value[i];
    }
    out << ']';
}

void dump_vector(std::ostream &out, const std::vector<bool> &value) {
    out << '[';
    for (size_t i = 0; i < value.size(); ++i) {
        if (i) out << ',';
        out << (value[i] ? "true" : "false");
    }
    out << ']';
}

void dump_vector(std::ostream &out, const std::vector<std::vector<double>> &value) {
    out << '[';
    for (size_t i = 0; i < value.size(); ++i) {
        if (i) out << ',';
        dump_vector(out, value[i]);
    }
    out << ']';
}

void dump_vector(std::ostream &out, const std::vector<std::vector<int>> &value) {
    out << '[';
    for (size_t i = 0; i < value.size(); ++i) {
        if (i) out << ',';
        dump_vector(out, value[i]);
    }
    out << ']';
}

void dump_vector(std::ostream &out, const std::vector<std::vector<std::vector<double>>> &value) {
    out << '[';
    for (size_t i = 0; i < value.size(); ++i) {
        if (i) out << ',';
        dump_vector(out, value[i]);
    }
    out << ']';
}

#define ADD_SCALAR(map, obj, field) map.emplace(#field, [&](const Json &v) { assign_scalar(obj.field, v); })
#define ADD_ARRAY(map, obj, field) map.emplace(#field, [&](const Json &v) { assign_vector(obj.field, v); })

void populate_ctx_setters(XFoilState &ctx, std::unordered_map<std::string, std::function<void(const Json &)>> &setters) {
    ADD_SCALAR(setters, ctx, LALFA);
    ADD_SCALAR(setters, ctx, LBLINI);
    ADD_SCALAR(setters, ctx, CL);
    ADD_SCALAR(setters, ctx, CM);
    ADD_SCALAR(setters, ctx, CD);
    ADD_SCALAR(setters, ctx, CDP);
    ADD_SCALAR(setters, ctx, CDF);
    ADD_SCALAR(setters, ctx, CL_ALF);
    ADD_SCALAR(setters, ctx, CL_MSQ);
    ADD_SCALAR(setters, ctx, CLSPEC);
    ADD_SCALAR(setters, ctx, MINF);
    ADD_SCALAR(setters, ctx, MINF1);
    ADD_SCALAR(setters, ctx, MINF_CL);
    ADD_SCALAR(setters, ctx, QINF);
    ADD_SCALAR(setters, ctx, TKLAM);
    ADD_SCALAR(setters, ctx, TKL_MSQ);
    ADD_SCALAR(setters, ctx, GAMMA);
    ADD_SCALAR(setters, ctx, GAMM1);
    ADD_SCALAR(setters, ctx, HVRAT);
    ADD_SCALAR(setters, ctx, REINF);
    ADD_SCALAR(setters, ctx, REINF1);
    ADD_SCALAR(setters, ctx, REINF_CL);
    ADD_SCALAR(setters, ctx, RETYP);
    ADD_SCALAR(setters, ctx, MATYP);
    ADD_SCALAR(setters, ctx, CPSTAR);
    ADD_SCALAR(setters, ctx, QSTAR);
    ADD_SCALAR(setters, ctx, IDAMP);
    ADD_SCALAR(setters, ctx, SLE);
    ADD_SCALAR(setters, ctx, XLE);
    ADD_SCALAR(setters, ctx, YLE);
    ADD_SCALAR(setters, ctx, XTE);
    ADD_SCALAR(setters, ctx, YTE);
    ADD_SCALAR(setters, ctx, SST);
    ADD_SCALAR(setters, ctx, SST_GO);
    ADD_SCALAR(setters, ctx, SST_GP);
    ADD_SCALAR(setters, ctx, ALFA);
    ADD_SCALAR(setters, ctx, ADEG);
    ADD_SCALAR(setters, ctx, DTOR);
    ADD_SCALAR(setters, ctx, AVISC);
    ADD_SCALAR(setters, ctx, MVISC);
    ADD_SCALAR(setters, ctx, XCMREF);
    ADD_SCALAR(setters, ctx, YCMREF);
    ADD_SCALAR(setters, ctx, N);
    ADD_SCALAR(setters, ctx, IST);
    ADD_SCALAR(setters, ctx, NW);
    ADD_SCALAR(setters, ctx, NPAN);
    ADD_SCALAR(setters, ctx, NB);
    ADD_SCALAR(setters, ctx, LCLOCK);

    ADD_ARRAY(setters, ctx, X);
    ADD_ARRAY(setters, ctx, Y);
    ADD_ARRAY(setters, ctx, XP);
    ADD_ARRAY(setters, ctx, YP);
    ADD_ARRAY(setters, ctx, S);
    ADD_ARRAY(setters, ctx, SNEW);
    ADD_ARRAY(setters, ctx, W1);
    ADD_ARRAY(setters, ctx, W2);
    ADD_ARRAY(setters, ctx, W3);
    ADD_ARRAY(setters, ctx, W4);
    ADD_ARRAY(setters, ctx, W5);
    ADD_ARRAY(setters, ctx, W6);
    ADD_ARRAY(setters, ctx, W7);
    ADD_ARRAY(setters, ctx, W8);
    ADD_ARRAY(setters, ctx, XB);
    ADD_ARRAY(setters, ctx, YB);
    ADD_ARRAY(setters, ctx, XBP);
    ADD_ARRAY(setters, ctx, YBP);
    ADD_ARRAY(setters, ctx, SB);

    ADD_SCALAR(setters, ctx, SBLE);
    ADD_SCALAR(setters, ctx, CHORDB);
    ADD_SCALAR(setters, ctx, AREAB);
    ADD_SCALAR(setters, ctx, RADBLE);
    ADD_SCALAR(setters, ctx, ANGBTE);
    ADD_SCALAR(setters, ctx, EI11BA);
    ADD_SCALAR(setters, ctx, EI22BA);
    ADD_SCALAR(setters, ctx, APX1BA);
    ADD_SCALAR(setters, ctx, APX2BA);
    ADD_SCALAR(setters, ctx, EI11BT);
    ADD_SCALAR(setters, ctx, EI22BT);
    ADD_SCALAR(setters, ctx, APX1BT);
    ADD_SCALAR(setters, ctx, APX2BT);
    ADD_SCALAR(setters, ctx, THICKB);
    ADD_SCALAR(setters, ctx, CAMBRB);

    ADD_ARRAY(setters, ctx, XSSI);
    ADD_ARRAY(setters, ctx, UEDG);
    ADD_ARRAY(setters, ctx, UINV);
    ADD_ARRAY(setters, ctx, UINV_A);
    ADD_ARRAY(setters, ctx, MASS);
    ADD_ARRAY(setters, ctx, THET);
    ADD_ARRAY(setters, ctx, DSTR);
    ADD_ARRAY(setters, ctx, CTAU);
    ADD_ARRAY(setters, ctx, DELT);
    ADD_ARRAY(setters, ctx, TSTR);
    ADD_ARRAY(setters, ctx, USLP);
    ADD_ARRAY(setters, ctx, GUXQ);
    ADD_ARRAY(setters, ctx, GUXD);
    ADD_ARRAY(setters, ctx, TAU);
    ADD_ARRAY(setters, ctx, DIS);
    ADD_ARRAY(setters, ctx, CTQ);
    ADD_ARRAY(setters, ctx, VTI);

    ADD_ARRAY(setters, ctx, ACRIT);
    ADD_ARRAY(setters, ctx, XSTRIP);
    ADD_ARRAY(setters, ctx, XOCTR);
    ADD_ARRAY(setters, ctx, YOCTR);
    ADD_ARRAY(setters, ctx, XSSITR);
    ADD_ARRAY(setters, ctx, TINDEX);
    ADD_ARRAY(setters, ctx, IBLTE);
    ADD_ARRAY(setters, ctx, NBL);
    ADD_ARRAY(setters, ctx, IPAN);
    ADD_ARRAY(setters, ctx, ISYS);
    ADD_SCALAR(setters, ctx, NSYS);
    ADD_ARRAY(setters, ctx, ITRAN);
    ADD_ARRAY(setters, ctx, TFORCE);

    ADD_ARRAY(setters, ctx, WGAP);
    ADD_SCALAR(setters, ctx, DWTE);
    ADD_SCALAR(setters, ctx, ANTE);
    ADD_SCALAR(setters, ctx, DSTE);
    ADD_SCALAR(setters, ctx, ASTE);
    ADD_SCALAR(setters, ctx, WAKLEN);
    ADD_SCALAR(setters, ctx, CHORD);
    ADD_SCALAR(setters, ctx, YIMAGE);
    ADD_SCALAR(setters, ctx, SHARP);

    ADD_ARRAY(setters, ctx, DIJ);
    ADD_ARRAY(setters, ctx, AIJ);
    ADD_ARRAY(setters, ctx, BIJ);
    ADD_ARRAY(setters, ctx, CIJ);
    ADD_ARRAY(setters, ctx, VM);
    ADD_ARRAY(setters, ctx, VA);
    ADD_ARRAY(setters, ctx, VB);
    ADD_ARRAY(setters, ctx, VDEL);
    ADD_ARRAY(setters, ctx, VZ);
    ADD_ARRAY(setters, ctx, AIJPIV);

    ADD_ARRAY(setters, ctx, QINV);
    ADD_ARRAY(setters, ctx, QVIS);
    ADD_ARRAY(setters, ctx, CPI);
    ADD_ARRAY(setters, ctx, CPV);
    ADD_ARRAY(setters, ctx, QINVU);
    ADD_ARRAY(setters, ctx, QINV_A);

    ADD_ARRAY(setters, ctx, GAM);
    ADD_ARRAY(setters, ctx, GAMU);
    ADD_ARRAY(setters, ctx, GAM_A);
    ADD_ARRAY(setters, ctx, SIG);
    ADD_ARRAY(setters, ctx, NX);
    ADD_ARRAY(setters, ctx, NY);
    ADD_ARRAY(setters, ctx, APANEL);
    ADD_SCALAR(setters, ctx, GAMTE);
    ADD_SCALAR(setters, ctx, GAMTE_A);
    ADD_SCALAR(setters, ctx, SIGTE);
    ADD_SCALAR(setters, ctx, SIGTE_A);

    ADD_ARRAY(setters, ctx, DZDG);
    ADD_ARRAY(setters, ctx, DZDN);
    ADD_ARRAY(setters, ctx, DZDM);
    ADD_ARRAY(setters, ctx, DQDG);
    ADD_ARRAY(setters, ctx, DQDM);
    ADD_SCALAR(setters, ctx, QTAN1);
    ADD_SCALAR(setters, ctx, QTAN2);
    ADD_SCALAR(setters, ctx, Z_QINF);
    ADD_SCALAR(setters, ctx, Z_ALFA);
    ADD_SCALAR(setters, ctx, Z_QDOF0);
    ADD_SCALAR(setters, ctx, Z_QDOF1);
    ADD_SCALAR(setters, ctx, Z_QDOF2);
    ADD_SCALAR(setters, ctx, Z_QDOF3);

    ADD_ARRAY(setters, ctx, QF0);
    ADD_ARRAY(setters, ctx, QF1);
    ADD_ARRAY(setters, ctx, QF2);
    ADD_ARRAY(setters, ctx, QF3);

    ADD_SCALAR(setters, ctx, PI);
    ADD_SCALAR(setters, ctx, HOPI);
    ADD_SCALAR(setters, ctx, QOPI);
    ADD_SCALAR(setters, ctx, RMSBL);
    ADD_SCALAR(setters, ctx, RMXBL);
    ADD_SCALAR(setters, ctx, RLX);
    ADD_SCALAR(setters, ctx, VACCEL);
    ADD_SCALAR(setters, ctx, IMXBL);
    ADD_SCALAR(setters, ctx, ISMXBL);
    ADD_SCALAR(setters, ctx, VMXBL);
    ADD_SCALAR(setters, ctx, LQAIJ);
    ADD_SCALAR(setters, ctx, LADIJ);
    ADD_SCALAR(setters, ctx, LWDIJ);
    ADD_SCALAR(setters, ctx, LWAKE);
    ADD_SCALAR(setters, ctx, LGAMU);
    ADD_SCALAR(setters, ctx, LVISC);
    ADD_SCALAR(setters, ctx, LVCONV);
    ADD_SCALAR(setters, ctx, LFLAP);
    ADD_SCALAR(setters, ctx, LIMAGE);
    ADD_SCALAR(setters, ctx, LQINU);
    ADD_SCALAR(setters, ctx, LQSPEC);
    ADD_SCALAR(setters, ctx, LGSAME);
    ADD_SCALAR(setters, ctx, LSCINI);
    ADD_SCALAR(setters, ctx, LIPAN);
    ADD_SCALAR(setters, ctx, AWAKE);
    ADD_SCALAR(setters, ctx, XOF);
    ADD_SCALAR(setters, ctx, YOF);
    ADD_SCALAR(setters, ctx, HMOM);
    ADD_SCALAR(setters, ctx, HFX);
    ADD_SCALAR(setters, ctx, HFY);
    ADD_SCALAR(setters, ctx, CVPAR);
    ADD_SCALAR(setters, ctx, CTERAT);
    ADD_SCALAR(setters, ctx, CTRRAT);
    ADD_SCALAR(setters, ctx, XSREF1);
    ADD_SCALAR(setters, ctx, XSREF2);
    ADD_SCALAR(setters, ctx, XPREF1);
    ADD_SCALAR(setters, ctx, XPREF2);
    ADD_SCALAR(setters, ctx, XBF);
    ADD_SCALAR(setters, ctx, YBF);
    ADD_SCALAR(setters, ctx, LBFLAP);
    ADD_SCALAR(setters, ctx, NAME);
    ADD_SCALAR(setters, ctx, NNAME);
    ADD_ARRAY(setters, ctx, HTARG);
}

void populate_bl_setters(XBlState &bl, std::unordered_map<std::string, std::function<void(const Json &)>> &setters) {
    ADD_ARRAY(setters, bl, COM1);
    ADD_ARRAY(setters, bl, COM2);
    ADD_ARRAY(setters, bl, C1SAV);
    ADD_ARRAY(setters, bl, C2SAV);

    ADD_SCALAR(setters, bl, SIMI);
    ADD_SCALAR(setters, bl, TRAN);
    ADD_SCALAR(setters, bl, TURB);
    ADD_SCALAR(setters, bl, WAKE);
    ADD_SCALAR(setters, bl, TRFORC);
    ADD_SCALAR(setters, bl, TRFREE);
    ADD_SCALAR(setters, bl, IDAMPV);

    ADD_ARRAY(setters, bl, VS1);
    ADD_ARRAY(setters, bl, VS2);
    ADD_ARRAY(setters, bl, VSREZ);
    ADD_ARRAY(setters, bl, VSR);
    ADD_ARRAY(setters, bl, VSM);
    ADD_ARRAY(setters, bl, VSX);

    ADD_SCALAR(setters, bl, SCCON);
    ADD_SCALAR(setters, bl, GACON);
    ADD_SCALAR(setters, bl, GBCON);
    ADD_SCALAR(setters, bl, GCCON);
    ADD_SCALAR(setters, bl, DLCON);
    ADD_SCALAR(setters, bl, CTRCON);
    ADD_SCALAR(setters, bl, CTRCEX);
    ADD_SCALAR(setters, bl, DUXCON);
    ADD_SCALAR(setters, bl, CTCON);
    ADD_SCALAR(setters, bl, CFFAC);

    ADD_SCALAR(setters, bl, X1);
    ADD_SCALAR(setters, bl, U1);
    ADD_SCALAR(setters, bl, T1);
    ADD_SCALAR(setters, bl, D1);
    ADD_SCALAR(setters, bl, S1);
    ADD_SCALAR(setters, bl, AMPL1);
    ADD_SCALAR(setters, bl, U1_UEI);
    ADD_SCALAR(setters, bl, U1_MS);
    ADD_SCALAR(setters, bl, DW1);
    ADD_SCALAR(setters, bl, H1);
    ADD_SCALAR(setters, bl, H1_T1);
    ADD_SCALAR(setters, bl, H1_D1);
    ADD_SCALAR(setters, bl, M1);
    ADD_SCALAR(setters, bl, M1_U1);
    ADD_SCALAR(setters, bl, M1_MS);
    ADD_SCALAR(setters, bl, R1);
    ADD_SCALAR(setters, bl, R1_U1);
    ADD_SCALAR(setters, bl, R1_MS);
    ADD_SCALAR(setters, bl, V1);
    ADD_SCALAR(setters, bl, V1_U1);
    ADD_SCALAR(setters, bl, V1_MS);
    ADD_SCALAR(setters, bl, V1_RE);
    ADD_SCALAR(setters, bl, HK1);
    ADD_SCALAR(setters, bl, HK1_U1);
    ADD_SCALAR(setters, bl, HK1_T1);
    ADD_SCALAR(setters, bl, HK1_D1);
    ADD_SCALAR(setters, bl, HK1_MS);
    ADD_SCALAR(setters, bl, HS1);
    ADD_SCALAR(setters, bl, HS1_U1);
    ADD_SCALAR(setters, bl, HS1_T1);
    ADD_SCALAR(setters, bl, HS1_D1);
    ADD_SCALAR(setters, bl, HS1_MS);
    ADD_SCALAR(setters, bl, HS1_RE);
    ADD_SCALAR(setters, bl, HC1);
    ADD_SCALAR(setters, bl, HC1_U1);
    ADD_SCALAR(setters, bl, HC1_T1);
    ADD_SCALAR(setters, bl, HC1_D1);
    ADD_SCALAR(setters, bl, HC1_MS);
    ADD_SCALAR(setters, bl, RT1);
    ADD_SCALAR(setters, bl, RT1_U1);
    ADD_SCALAR(setters, bl, RT1_T1);
    ADD_SCALAR(setters, bl, RT1_MS);
    ADD_SCALAR(setters, bl, RT1_RE);
    ADD_SCALAR(setters, bl, CF1);
    ADD_SCALAR(setters, bl, CF1_U1);
    ADD_SCALAR(setters, bl, CF1_T1);
    ADD_SCALAR(setters, bl, CF1_D1);
    ADD_SCALAR(setters, bl, CF1_MS);
    ADD_SCALAR(setters, bl, CF1_RE);
    ADD_SCALAR(setters, bl, DI1);
    ADD_SCALAR(setters, bl, DI1_U1);
    ADD_SCALAR(setters, bl, DI1_T1);
    ADD_SCALAR(setters, bl, DI1_D1);
    ADD_SCALAR(setters, bl, DI1_S1);
    ADD_SCALAR(setters, bl, DI1_MS);
    ADD_SCALAR(setters, bl, DI1_RE);
    ADD_SCALAR(setters, bl, US1);
    ADD_SCALAR(setters, bl, US1_U1);
    ADD_SCALAR(setters, bl, US1_T1);
    ADD_SCALAR(setters, bl, US1_D1);
    ADD_SCALAR(setters, bl, US1_MS);
    ADD_SCALAR(setters, bl, US1_RE);
    ADD_SCALAR(setters, bl, CQ1);
    ADD_SCALAR(setters, bl, CQ1_U1);
    ADD_SCALAR(setters, bl, CQ1_T1);
    ADD_SCALAR(setters, bl, CQ1_D1);
    ADD_SCALAR(setters, bl, CQ1_MS);
    ADD_SCALAR(setters, bl, CQ1_RE);
    ADD_SCALAR(setters, bl, DE1);
    ADD_SCALAR(setters, bl, DE1_U1);
    ADD_SCALAR(setters, bl, DE1_T1);
    ADD_SCALAR(setters, bl, DE1_D1);
    ADD_SCALAR(setters, bl, DE1_MS);

    ADD_SCALAR(setters, bl, X2);
    ADD_SCALAR(setters, bl, U2);
    ADD_SCALAR(setters, bl, T2);
    ADD_SCALAR(setters, bl, D2);
    ADD_SCALAR(setters, bl, S2);
    ADD_SCALAR(setters, bl, AMPL2);
    ADD_SCALAR(setters, bl, U2_UEI);
    ADD_SCALAR(setters, bl, U2_MS);
    ADD_SCALAR(setters, bl, DW2);
    ADD_SCALAR(setters, bl, H2);
    ADD_SCALAR(setters, bl, H2_T2);
    ADD_SCALAR(setters, bl, H2_D2);
    ADD_SCALAR(setters, bl, M2);
    ADD_SCALAR(setters, bl, M2_U2);
    ADD_SCALAR(setters, bl, M2_MS);
    ADD_SCALAR(setters, bl, R2);
    ADD_SCALAR(setters, bl, R2_U2);
    ADD_SCALAR(setters, bl, R2_MS);
    ADD_SCALAR(setters, bl, V2);
    ADD_SCALAR(setters, bl, V2_U2);
    ADD_SCALAR(setters, bl, V2_MS);
    ADD_SCALAR(setters, bl, V2_RE);
    ADD_SCALAR(setters, bl, HK2);
    ADD_SCALAR(setters, bl, HK2_U2);
    ADD_SCALAR(setters, bl, HK2_T2);
    ADD_SCALAR(setters, bl, HK2_D2);
    ADD_SCALAR(setters, bl, HK2_MS);
    ADD_SCALAR(setters, bl, HS2);
    ADD_SCALAR(setters, bl, HS2_U2);
    ADD_SCALAR(setters, bl, HS2_T2);
    ADD_SCALAR(setters, bl, HS2_D2);
    ADD_SCALAR(setters, bl, HS2_MS);
    ADD_SCALAR(setters, bl, HS2_RE);
    ADD_SCALAR(setters, bl, HC2);
    ADD_SCALAR(setters, bl, HC2_U2);
    ADD_SCALAR(setters, bl, HC2_T2);
    ADD_SCALAR(setters, bl, HC2_D2);
    ADD_SCALAR(setters, bl, HC2_MS);
    ADD_SCALAR(setters, bl, RT2);
    ADD_SCALAR(setters, bl, RT2_U2);
    ADD_SCALAR(setters, bl, RT2_T2);
    ADD_SCALAR(setters, bl, RT2_MS);
    ADD_SCALAR(setters, bl, RT2_RE);
    ADD_SCALAR(setters, bl, CF2);
    ADD_SCALAR(setters, bl, CF2_HK2);
    ADD_SCALAR(setters, bl, CF2_M2);
    ADD_SCALAR(setters, bl, CF2_RT2);
    ADD_SCALAR(setters, bl, CF2_U2);
    ADD_SCALAR(setters, bl, CF2_T2);
    ADD_SCALAR(setters, bl, CF2_D2);
    ADD_SCALAR(setters, bl, CF2_MS);
    ADD_SCALAR(setters, bl, CF2_RE);
    ADD_SCALAR(setters, bl, DI2);
    ADD_SCALAR(setters, bl, DI2_U2);
    ADD_SCALAR(setters, bl, DI2_T2);
    ADD_SCALAR(setters, bl, DI2_D2);
    ADD_SCALAR(setters, bl, DI2_S2);
    ADD_SCALAR(setters, bl, DI2_MS);
    ADD_SCALAR(setters, bl, DI2_RE);
    ADD_SCALAR(setters, bl, US2);
    ADD_SCALAR(setters, bl, US2_U2);
    ADD_SCALAR(setters, bl, US2_T2);
    ADD_SCALAR(setters, bl, US2_D2);
    ADD_SCALAR(setters, bl, US2_MS);
    ADD_SCALAR(setters, bl, US2_RE);
    ADD_SCALAR(setters, bl, CQ2);
    ADD_SCALAR(setters, bl, CQ2_U2);
    ADD_SCALAR(setters, bl, CQ2_T2);
    ADD_SCALAR(setters, bl, CQ2_D2);
    ADD_SCALAR(setters, bl, CQ2_MS);
    ADD_SCALAR(setters, bl, CQ2_RE);
    ADD_SCALAR(setters, bl, DE2);
    ADD_SCALAR(setters, bl, DE2_U2);
    ADD_SCALAR(setters, bl, DE2_T2);
    ADD_SCALAR(setters, bl, DE2_D2);
    ADD_SCALAR(setters, bl, DE2_MS);

    ADD_SCALAR(setters, bl, CFM);
    ADD_SCALAR(setters, bl, CFM_HKA);
    ADD_SCALAR(setters, bl, CFM_MA);
    ADD_SCALAR(setters, bl, CFM_MS);
    ADD_SCALAR(setters, bl, CFM_RE);
    ADD_SCALAR(setters, bl, CFM_RTA);
    ADD_SCALAR(setters, bl, CFM_U1);
    ADD_SCALAR(setters, bl, CFM_T1);
    ADD_SCALAR(setters, bl, CFM_D1);
    ADD_SCALAR(setters, bl, CFM_U2);
    ADD_SCALAR(setters, bl, CFM_T2);
    ADD_SCALAR(setters, bl, CFM_D2);
    ADD_SCALAR(setters, bl, XT);
    ADD_SCALAR(setters, bl, XT_A1);
    ADD_SCALAR(setters, bl, XT_A2);
    ADD_SCALAR(setters, bl, XT_MS);
    ADD_SCALAR(setters, bl, XT_RE);
    ADD_SCALAR(setters, bl, XT_XF);
    ADD_SCALAR(setters, bl, XT_X1);
    ADD_SCALAR(setters, bl, XT_T1);
    ADD_SCALAR(setters, bl, XT_D1);
    ADD_SCALAR(setters, bl, XT_U1);
    ADD_SCALAR(setters, bl, XT_X2);
    ADD_SCALAR(setters, bl, XT_T2);
    ADD_SCALAR(setters, bl, XT_D2);
    ADD_SCALAR(setters, bl, XT_U2);
    ADD_SCALAR(setters, bl, DWTE);
    ADD_SCALAR(setters, bl, QINFBL);
    ADD_SCALAR(setters, bl, TKBL);
    ADD_SCALAR(setters, bl, TKBL_MS);
    ADD_SCALAR(setters, bl, RSTBL);
    ADD_SCALAR(setters, bl, RSTBL_MS);
    ADD_SCALAR(setters, bl, HSTINV);
    ADD_SCALAR(setters, bl, HSTINV_MS);
    ADD_SCALAR(setters, bl, REYBL);
    ADD_SCALAR(setters, bl, REYBL_MS);
    ADD_SCALAR(setters, bl, REYBL_RE);
    ADD_SCALAR(setters, bl, GAMBL);
    ADD_SCALAR(setters, bl, GM1BL);
    ADD_SCALAR(setters, bl, HVRAT);
    ADD_SCALAR(setters, bl, BULE);
    ADD_SCALAR(setters, bl, XIFORC);
    ADD_SCALAR(setters, bl, AMCRIT);
}

#undef ADD_SCALAR
#undef ADD_ARRAY

void dump_field(std::ostream &out, const std::string &name, const std::string &value, bool &first) {
    if (!first) {
        out << ',';
    }
    first = false;
    write_json_string(out, name);
    out << ':';
    write_json_string(out, value);
}

template <typename T>
void dump_field(std::ostream &out, const std::string &name, const T &value, bool &first) {
    if (!first) {
        out << ',';
    }
    first = false;
    write_json_string(out, name);
    out << ':';
    out << std::setprecision(17) << value;
}

void dump_field(std::ostream &out, const std::string &name, bool value, bool &first) {
    if (!first) {
        out << ',';
    }
    first = false;
    write_json_string(out, name);
    out << ':' << (value ? "true" : "false");
}

template <typename Vec>
void dump_field_vec(std::ostream &out, const std::string &name, const Vec &value, bool &first) {
    if (!first) {
        out << ',';
    }
    first = false;
    write_json_string(out, name);
    out << ':';
    dump_vector(out, value);
}

void dump_ctx(std::ostream &out, const XFoilState &ctx) {
    bool first = true;
    out << '{';
    dump_field(out, "LALFA", ctx.LALFA, first);
    dump_field(out, "LBLINI", ctx.LBLINI, first);
    dump_field(out, "CL", ctx.CL, first);
    dump_field(out, "CM", ctx.CM, first);
    dump_field(out, "CD", ctx.CD, first);
    dump_field(out, "CDP", ctx.CDP, first);
    dump_field(out, "CDF", ctx.CDF, first);
    dump_field(out, "CL_ALF", ctx.CL_ALF, first);
    dump_field(out, "CL_MSQ", ctx.CL_MSQ, first);
    dump_field(out, "CLSPEC", ctx.CLSPEC, first);
    dump_field(out, "MINF", ctx.MINF, first);
    dump_field(out, "MINF1", ctx.MINF1, first);
    dump_field(out, "MINF_CL", ctx.MINF_CL, first);
    dump_field(out, "QINF", ctx.QINF, first);
    dump_field(out, "TKLAM", ctx.TKLAM, first);
    dump_field(out, "TKL_MSQ", ctx.TKL_MSQ, first);
    dump_field(out, "GAMMA", ctx.GAMMA, first);
    dump_field(out, "GAMM1", ctx.GAMM1, first);
    dump_field(out, "HVRAT", ctx.HVRAT, first);
    dump_field(out, "REINF", ctx.REINF, first);
    dump_field(out, "REINF1", ctx.REINF1, first);
    dump_field(out, "REINF_CL", ctx.REINF_CL, first);
    dump_field(out, "RETYP", ctx.RETYP, first);
    dump_field(out, "MATYP", ctx.MATYP, first);
    dump_field(out, "CPSTAR", ctx.CPSTAR, first);
    dump_field(out, "QSTAR", ctx.QSTAR, first);
    dump_field(out, "IDAMP", ctx.IDAMP, first);
    dump_field(out, "SLE", ctx.SLE, first);
    dump_field(out, "XLE", ctx.XLE, first);
    dump_field(out, "YLE", ctx.YLE, first);
    dump_field(out, "XTE", ctx.XTE, first);
    dump_field(out, "YTE", ctx.YTE, first);
    dump_field(out, "SST", ctx.SST, first);
    dump_field(out, "SST_GO", ctx.SST_GO, first);
    dump_field(out, "SST_GP", ctx.SST_GP, first);
    dump_field(out, "ALFA", ctx.ALFA, first);
    dump_field(out, "ADEG", ctx.ADEG, first);
    dump_field(out, "DTOR", ctx.DTOR, first);
    dump_field(out, "AVISC", ctx.AVISC, first);
    dump_field(out, "MVISC", ctx.MVISC, first);
    dump_field(out, "XCMREF", ctx.XCMREF, first);
    dump_field(out, "YCMREF", ctx.YCMREF, first);
    dump_field(out, "N", ctx.N, first);
    dump_field(out, "IST", ctx.IST, first);
    dump_field(out, "NW", ctx.NW, first);
    dump_field(out, "NPAN", ctx.NPAN, first);
    dump_field(out, "NB", ctx.NB, first);
    dump_field(out, "LCLOCK", ctx.LCLOCK, first);

    dump_field_vec(out, "X", ctx.X, first);
    dump_field_vec(out, "Y", ctx.Y, first);
    dump_field_vec(out, "XP", ctx.XP, first);
    dump_field_vec(out, "YP", ctx.YP, first);
    dump_field_vec(out, "S", ctx.S, first);
    dump_field_vec(out, "SNEW", ctx.SNEW, first);
    dump_field_vec(out, "W1", ctx.W1, first);
    dump_field_vec(out, "W2", ctx.W2, first);
    dump_field_vec(out, "W3", ctx.W3, first);
    dump_field_vec(out, "W4", ctx.W4, first);
    dump_field_vec(out, "W5", ctx.W5, first);
    dump_field_vec(out, "W6", ctx.W6, first);
    dump_field_vec(out, "W7", ctx.W7, first);
    dump_field_vec(out, "W8", ctx.W8, first);
    dump_field_vec(out, "XB", ctx.XB, first);
    dump_field_vec(out, "YB", ctx.YB, first);
    dump_field_vec(out, "XBP", ctx.XBP, first);
    dump_field_vec(out, "YBP", ctx.YBP, first);
    dump_field_vec(out, "SB", ctx.SB, first);

    dump_field(out, "SBLE", ctx.SBLE, first);
    dump_field(out, "CHORDB", ctx.CHORDB, first);
    dump_field(out, "AREAB", ctx.AREAB, first);
    dump_field(out, "RADBLE", ctx.RADBLE, first);
    dump_field(out, "ANGBTE", ctx.ANGBTE, first);
    dump_field(out, "EI11BA", ctx.EI11BA, first);
    dump_field(out, "EI22BA", ctx.EI22BA, first);
    dump_field(out, "APX1BA", ctx.APX1BA, first);
    dump_field(out, "APX2BA", ctx.APX2BA, first);
    dump_field(out, "EI11BT", ctx.EI11BT, first);
    dump_field(out, "EI22BT", ctx.EI22BT, first);
    dump_field(out, "APX1BT", ctx.APX1BT, first);
    dump_field(out, "APX2BT", ctx.APX2BT, first);
    dump_field(out, "THICKB", ctx.THICKB, first);
    dump_field(out, "CAMBRB", ctx.CAMBRB, first);

    dump_field_vec(out, "XSSI", ctx.XSSI, first);
    dump_field_vec(out, "UEDG", ctx.UEDG, first);
    dump_field_vec(out, "UINV", ctx.UINV, first);
    dump_field_vec(out, "UINV_A", ctx.UINV_A, first);
    dump_field_vec(out, "MASS", ctx.MASS, first);
    dump_field_vec(out, "THET", ctx.THET, first);
    dump_field_vec(out, "DSTR", ctx.DSTR, first);
    dump_field_vec(out, "CTAU", ctx.CTAU, first);
    dump_field_vec(out, "DELT", ctx.DELT, first);
    dump_field_vec(out, "TSTR", ctx.TSTR, first);
    dump_field_vec(out, "USLP", ctx.USLP, first);
    dump_field_vec(out, "GUXQ", ctx.GUXQ, first);
    dump_field_vec(out, "GUXD", ctx.GUXD, first);
    dump_field_vec(out, "TAU", ctx.TAU, first);
    dump_field_vec(out, "DIS", ctx.DIS, first);
    dump_field_vec(out, "CTQ", ctx.CTQ, first);
    dump_field_vec(out, "VTI", ctx.VTI, first);

    dump_field_vec(out, "ACRIT", ctx.ACRIT, first);
    dump_field_vec(out, "XSTRIP", ctx.XSTRIP, first);
    dump_field_vec(out, "XOCTR", ctx.XOCTR, first);
    dump_field_vec(out, "YOCTR", ctx.YOCTR, first);
    dump_field_vec(out, "XSSITR", ctx.XSSITR, first);
    dump_field_vec(out, "TINDEX", ctx.TINDEX, first);
    dump_field_vec(out, "IBLTE", ctx.IBLTE, first);
    dump_field_vec(out, "NBL", ctx.NBL, first);
    dump_field_vec(out, "IPAN", ctx.IPAN, first);
    dump_field_vec(out, "ISYS", ctx.ISYS, first);
    dump_field(out, "NSYS", ctx.NSYS, first);
    dump_field_vec(out, "ITRAN", ctx.ITRAN, first);
    dump_field_vec(out, "TFORCE", ctx.TFORCE, first);

    dump_field_vec(out, "WGAP", ctx.WGAP, first);
    dump_field(out, "DWTE", ctx.DWTE, first);
    dump_field(out, "ANTE", ctx.ANTE, first);
    dump_field(out, "DSTE", ctx.DSTE, first);
    dump_field(out, "ASTE", ctx.ASTE, first);
    dump_field(out, "WAKLEN", ctx.WAKLEN, first);
    dump_field(out, "CHORD", ctx.CHORD, first);
    dump_field(out, "YIMAGE", ctx.YIMAGE, first);
    dump_field(out, "SHARP", ctx.SHARP, first);

    dump_field_vec(out, "DIJ", ctx.DIJ, first);
    dump_field_vec(out, "AIJ", ctx.AIJ, first);
    dump_field_vec(out, "BIJ", ctx.BIJ, first);
    dump_field_vec(out, "CIJ", ctx.CIJ, first);
    dump_field_vec(out, "VM", ctx.VM, first);
    dump_field_vec(out, "VA", ctx.VA, first);
    dump_field_vec(out, "VB", ctx.VB, first);
    dump_field_vec(out, "VDEL", ctx.VDEL, first);
    dump_field_vec(out, "VZ", ctx.VZ, first);
    dump_field_vec(out, "AIJPIV", ctx.AIJPIV, first);

    dump_field_vec(out, "QINV", ctx.QINV, first);
    dump_field_vec(out, "QVIS", ctx.QVIS, first);
    dump_field_vec(out, "CPI", ctx.CPI, first);
    dump_field_vec(out, "CPV", ctx.CPV, first);
    dump_field_vec(out, "QINVU", ctx.QINVU, first);
    dump_field_vec(out, "QINV_A", ctx.QINV_A, first);

    dump_field_vec(out, "GAM", ctx.GAM, first);
    dump_field_vec(out, "GAMU", ctx.GAMU, first);
    dump_field_vec(out, "GAM_A", ctx.GAM_A, first);
    dump_field_vec(out, "SIG", ctx.SIG, first);
    dump_field_vec(out, "NX", ctx.NX, first);
    dump_field_vec(out, "NY", ctx.NY, first);
    dump_field_vec(out, "APANEL", ctx.APANEL, first);
    dump_field(out, "GAMTE", ctx.GAMTE, first);
    dump_field(out, "GAMTE_A", ctx.GAMTE_A, first);
    dump_field(out, "SIGTE", ctx.SIGTE, first);
    dump_field(out, "SIGTE_A", ctx.SIGTE_A, first);

    dump_field_vec(out, "DZDG", ctx.DZDG, first);
    dump_field_vec(out, "DZDN", ctx.DZDN, first);
    dump_field_vec(out, "DZDM", ctx.DZDM, first);
    dump_field_vec(out, "DQDG", ctx.DQDG, first);
    dump_field_vec(out, "DQDM", ctx.DQDM, first);
    dump_field(out, "QTAN1", ctx.QTAN1, first);
    dump_field(out, "QTAN2", ctx.QTAN2, first);
    dump_field(out, "Z_QINF", ctx.Z_QINF, first);
    dump_field(out, "Z_ALFA", ctx.Z_ALFA, first);
    dump_field(out, "Z_QDOF0", ctx.Z_QDOF0, first);
    dump_field(out, "Z_QDOF1", ctx.Z_QDOF1, first);
    dump_field(out, "Z_QDOF2", ctx.Z_QDOF2, first);
    dump_field(out, "Z_QDOF3", ctx.Z_QDOF3, first);

    dump_field_vec(out, "QF0", ctx.QF0, first);
    dump_field_vec(out, "QF1", ctx.QF1, first);
    dump_field_vec(out, "QF2", ctx.QF2, first);
    dump_field_vec(out, "QF3", ctx.QF3, first);

    dump_field(out, "PI", ctx.PI, first);
    dump_field(out, "HOPI", ctx.HOPI, first);
    dump_field(out, "QOPI", ctx.QOPI, first);
    dump_field(out, "RMSBL", ctx.RMSBL, first);
    dump_field(out, "RMXBL", ctx.RMXBL, first);
    dump_field(out, "RLX", ctx.RLX, first);
    dump_field(out, "VACCEL", ctx.VACCEL, first);
    dump_field(out, "IMXBL", ctx.IMXBL, first);
    dump_field(out, "ISMXBL", ctx.ISMXBL, first);
    dump_field(out, "VMXBL", ctx.VMXBL, first);
    dump_field(out, "LQAIJ", ctx.LQAIJ, first);
    dump_field(out, "LADIJ", ctx.LADIJ, first);
    dump_field(out, "LWDIJ", ctx.LWDIJ, first);
    dump_field(out, "LWAKE", ctx.LWAKE, first);
    dump_field(out, "LGAMU", ctx.LGAMU, first);
    dump_field(out, "LVISC", ctx.LVISC, first);
    dump_field(out, "LVCONV", ctx.LVCONV, first);
    dump_field(out, "LFLAP", ctx.LFLAP, first);
    dump_field(out, "LIMAGE", ctx.LIMAGE, first);
    dump_field(out, "LQINU", ctx.LQINU, first);
    dump_field(out, "LQSPEC", ctx.LQSPEC, first);
    dump_field(out, "LGSAME", ctx.LGSAME, first);
    dump_field(out, "LSCINI", ctx.LSCINI, first);
    dump_field(out, "LIPAN", ctx.LIPAN, first);
    dump_field(out, "AWAKE", ctx.AWAKE, first);
    dump_field(out, "XOF", ctx.XOF, first);
    dump_field(out, "YOF", ctx.YOF, first);
    dump_field(out, "HMOM", ctx.HMOM, first);
    dump_field(out, "HFX", ctx.HFX, first);
    dump_field(out, "HFY", ctx.HFY, first);
    dump_field(out, "CVPAR", ctx.CVPAR, first);
    dump_field(out, "CTERAT", ctx.CTERAT, first);
    dump_field(out, "CTRRAT", ctx.CTRRAT, first);
    dump_field(out, "XSREF1", ctx.XSREF1, first);
    dump_field(out, "XSREF2", ctx.XSREF2, first);
    dump_field(out, "XPREF1", ctx.XPREF1, first);
    dump_field(out, "XPREF2", ctx.XPREF2, first);
    dump_field(out, "XBF", ctx.XBF, first);
    dump_field(out, "YBF", ctx.YBF, first);
    dump_field(out, "LBFLAP", ctx.LBFLAP, first);
    dump_field(out, "NAME", ctx.NAME, first);
    dump_field(out, "NNAME", ctx.NNAME, first);
    dump_field_vec(out, "HTARG", ctx.HTARG, first);
    out << '}';
}

void dump_bl(std::ostream &out, const XBlState &bl) {
    bool first = true;
    out << '{';
    dump_field_vec(out, "COM1", bl.COM1, first);
    dump_field_vec(out, "COM2", bl.COM2, first);
    dump_field_vec(out, "C1SAV", bl.C1SAV, first);
    dump_field_vec(out, "C2SAV", bl.C2SAV, first);

    dump_field(out, "SIMI", bl.SIMI, first);
    dump_field(out, "TRAN", bl.TRAN, first);
    dump_field(out, "TURB", bl.TURB, first);
    dump_field(out, "WAKE", bl.WAKE, first);
    dump_field(out, "TRFORC", bl.TRFORC, first);
    dump_field(out, "TRFREE", bl.TRFREE, first);
    dump_field(out, "IDAMPV", bl.IDAMPV, first);

    dump_field_vec(out, "VS1", bl.VS1, first);
    dump_field_vec(out, "VS2", bl.VS2, first);
    dump_field_vec(out, "VSREZ", bl.VSREZ, first);
    dump_field_vec(out, "VSR", bl.VSR, first);
    dump_field_vec(out, "VSM", bl.VSM, first);
    dump_field_vec(out, "VSX", bl.VSX, first);

    dump_field(out, "SCCON", bl.SCCON, first);
    dump_field(out, "GACON", bl.GACON, first);
    dump_field(out, "GBCON", bl.GBCON, first);
    dump_field(out, "GCCON", bl.GCCON, first);
    dump_field(out, "DLCON", bl.DLCON, first);
    dump_field(out, "CTRCON", bl.CTRCON, first);
    dump_field(out, "CTRCEX", bl.CTRCEX, first);
    dump_field(out, "DUXCON", bl.DUXCON, first);
    dump_field(out, "CTCON", bl.CTCON, first);
    dump_field(out, "CFFAC", bl.CFFAC, first);

    dump_field(out, "X1", bl.X1, first);
    dump_field(out, "U1", bl.U1, first);
    dump_field(out, "T1", bl.T1, first);
    dump_field(out, "D1", bl.D1, first);
    dump_field(out, "S1", bl.S1, first);
    dump_field(out, "AMPL1", bl.AMPL1, first);
    dump_field(out, "U1_UEI", bl.U1_UEI, first);
    dump_field(out, "U1_MS", bl.U1_MS, first);
    dump_field(out, "DW1", bl.DW1, first);
    dump_field(out, "H1", bl.H1, first);
    dump_field(out, "H1_T1", bl.H1_T1, first);
    dump_field(out, "H1_D1", bl.H1_D1, first);
    dump_field(out, "M1", bl.M1, first);
    dump_field(out, "M1_U1", bl.M1_U1, first);
    dump_field(out, "M1_MS", bl.M1_MS, first);
    dump_field(out, "R1", bl.R1, first);
    dump_field(out, "R1_U1", bl.R1_U1, first);
    dump_field(out, "R1_MS", bl.R1_MS, first);
    dump_field(out, "V1", bl.V1, first);
    dump_field(out, "V1_U1", bl.V1_U1, first);
    dump_field(out, "V1_MS", bl.V1_MS, first);
    dump_field(out, "V1_RE", bl.V1_RE, first);
    dump_field(out, "HK1", bl.HK1, first);
    dump_field(out, "HK1_U1", bl.HK1_U1, first);
    dump_field(out, "HK1_T1", bl.HK1_T1, first);
    dump_field(out, "HK1_D1", bl.HK1_D1, first);
    dump_field(out, "HK1_MS", bl.HK1_MS, first);
    dump_field(out, "HS1", bl.HS1, first);
    dump_field(out, "HS1_U1", bl.HS1_U1, first);
    dump_field(out, "HS1_T1", bl.HS1_T1, first);
    dump_field(out, "HS1_D1", bl.HS1_D1, first);
    dump_field(out, "HS1_MS", bl.HS1_MS, first);
    dump_field(out, "HS1_RE", bl.HS1_RE, first);
    dump_field(out, "HC1", bl.HC1, first);
    dump_field(out, "HC1_U1", bl.HC1_U1, first);
    dump_field(out, "HC1_T1", bl.HC1_T1, first);
    dump_field(out, "HC1_D1", bl.HC1_D1, first);
    dump_field(out, "HC1_MS", bl.HC1_MS, first);
    dump_field(out, "RT1", bl.RT1, first);
    dump_field(out, "RT1_U1", bl.RT1_U1, first);
    dump_field(out, "RT1_T1", bl.RT1_T1, first);
    dump_field(out, "RT1_MS", bl.RT1_MS, first);
    dump_field(out, "RT1_RE", bl.RT1_RE, first);
    dump_field(out, "CF1", bl.CF1, first);
    dump_field(out, "CF1_U1", bl.CF1_U1, first);
    dump_field(out, "CF1_T1", bl.CF1_T1, first);
    dump_field(out, "CF1_D1", bl.CF1_D1, first);
    dump_field(out, "CF1_MS", bl.CF1_MS, first);
    dump_field(out, "CF1_RE", bl.CF1_RE, first);
    dump_field(out, "DI1", bl.DI1, first);
    dump_field(out, "DI1_U1", bl.DI1_U1, first);
    dump_field(out, "DI1_T1", bl.DI1_T1, first);
    dump_field(out, "DI1_D1", bl.DI1_D1, first);
    dump_field(out, "DI1_S1", bl.DI1_S1, first);
    dump_field(out, "DI1_MS", bl.DI1_MS, first);
    dump_field(out, "DI1_RE", bl.DI1_RE, first);
    dump_field(out, "US1", bl.US1, first);
    dump_field(out, "US1_U1", bl.US1_U1, first);
    dump_field(out, "US1_T1", bl.US1_T1, first);
    dump_field(out, "US1_D1", bl.US1_D1, first);
    dump_field(out, "US1_MS", bl.US1_MS, first);
    dump_field(out, "US1_RE", bl.US1_RE, first);
    dump_field(out, "CQ1", bl.CQ1, first);
    dump_field(out, "CQ1_U1", bl.CQ1_U1, first);
    dump_field(out, "CQ1_T1", bl.CQ1_T1, first);
    dump_field(out, "CQ1_D1", bl.CQ1_D1, first);
    dump_field(out, "CQ1_MS", bl.CQ1_MS, first);
    dump_field(out, "CQ1_RE", bl.CQ1_RE, first);
    dump_field(out, "DE1", bl.DE1, first);
    dump_field(out, "DE1_U1", bl.DE1_U1, first);
    dump_field(out, "DE1_T1", bl.DE1_T1, first);
    dump_field(out, "DE1_D1", bl.DE1_D1, first);
    dump_field(out, "DE1_MS", bl.DE1_MS, first);

    dump_field(out, "X2", bl.X2, first);
    dump_field(out, "U2", bl.U2, first);
    dump_field(out, "T2", bl.T2, first);
    dump_field(out, "D2", bl.D2, first);
    dump_field(out, "S2", bl.S2, first);
    dump_field(out, "AMPL2", bl.AMPL2, first);
    dump_field(out, "U2_UEI", bl.U2_UEI, first);
    dump_field(out, "U2_MS", bl.U2_MS, first);
    dump_field(out, "DW2", bl.DW2, first);
    dump_field(out, "H2", bl.H2, first);
    dump_field(out, "H2_T2", bl.H2_T2, first);
    dump_field(out, "H2_D2", bl.H2_D2, first);
    dump_field(out, "M2", bl.M2, first);
    dump_field(out, "M2_U2", bl.M2_U2, first);
    dump_field(out, "M2_MS", bl.M2_MS, first);
    dump_field(out, "R2", bl.R2, first);
    dump_field(out, "R2_U2", bl.R2_U2, first);
    dump_field(out, "R2_MS", bl.R2_MS, first);
    dump_field(out, "V2", bl.V2, first);
    dump_field(out, "V2_U2", bl.V2_U2, first);
    dump_field(out, "V2_MS", bl.V2_MS, first);
    dump_field(out, "V2_RE", bl.V2_RE, first);
    dump_field(out, "HK2", bl.HK2, first);
    dump_field(out, "HK2_U2", bl.HK2_U2, first);
    dump_field(out, "HK2_T2", bl.HK2_T2, first);
    dump_field(out, "HK2_D2", bl.HK2_D2, first);
    dump_field(out, "HK2_MS", bl.HK2_MS, first);
    dump_field(out, "HS2", bl.HS2, first);
    dump_field(out, "HS2_U2", bl.HS2_U2, first);
    dump_field(out, "HS2_T2", bl.HS2_T2, first);
    dump_field(out, "HS2_D2", bl.HS2_D2, first);
    dump_field(out, "HS2_MS", bl.HS2_MS, first);
    dump_field(out, "HS2_RE", bl.HS2_RE, first);
    dump_field(out, "HC2", bl.HC2, first);
    dump_field(out, "HC2_U2", bl.HC2_U2, first);
    dump_field(out, "HC2_T2", bl.HC2_T2, first);
    dump_field(out, "HC2_D2", bl.HC2_D2, first);
    dump_field(out, "HC2_MS", bl.HC2_MS, first);
    dump_field(out, "RT2", bl.RT2, first);
    dump_field(out, "RT2_U2", bl.RT2_U2, first);
    dump_field(out, "RT2_T2", bl.RT2_T2, first);
    dump_field(out, "RT2_MS", bl.RT2_MS, first);
    dump_field(out, "RT2_RE", bl.RT2_RE, first);
    dump_field(out, "CF2", bl.CF2, first);
    dump_field(out, "CF2_HK2", bl.CF2_HK2, first);
    dump_field(out, "CF2_M2", bl.CF2_M2, first);
    dump_field(out, "CF2_RT2", bl.CF2_RT2, first);
    dump_field(out, "CF2_U2", bl.CF2_U2, first);
    dump_field(out, "CF2_T2", bl.CF2_T2, first);
    dump_field(out, "CF2_D2", bl.CF2_D2, first);
    dump_field(out, "CF2_MS", bl.CF2_MS, first);
    dump_field(out, "CF2_RE", bl.CF2_RE, first);
    dump_field(out, "DI2", bl.DI2, first);
    dump_field(out, "DI2_U2", bl.DI2_U2, first);
    dump_field(out, "DI2_T2", bl.DI2_T2, first);
    dump_field(out, "DI2_D2", bl.DI2_D2, first);
    dump_field(out, "DI2_S2", bl.DI2_S2, first);
    dump_field(out, "DI2_MS", bl.DI2_MS, first);
    dump_field(out, "DI2_RE", bl.DI2_RE, first);
    dump_field(out, "US2", bl.US2, first);
    dump_field(out, "US2_U2", bl.US2_U2, first);
    dump_field(out, "US2_T2", bl.US2_T2, first);
    dump_field(out, "US2_D2", bl.US2_D2, first);
    dump_field(out, "US2_MS", bl.US2_MS, first);
    dump_field(out, "US2_RE", bl.US2_RE, first);
    dump_field(out, "CQ2", bl.CQ2, first);
    dump_field(out, "CQ2_U2", bl.CQ2_U2, first);
    dump_field(out, "CQ2_T2", bl.CQ2_T2, first);
    dump_field(out, "CQ2_D2", bl.CQ2_D2, first);
    dump_field(out, "CQ2_MS", bl.CQ2_MS, first);
    dump_field(out, "CQ2_RE", bl.CQ2_RE, first);
    dump_field(out, "DE2", bl.DE2, first);
    dump_field(out, "DE2_U2", bl.DE2_U2, first);
    dump_field(out, "DE2_T2", bl.DE2_T2, first);
    dump_field(out, "DE2_D2", bl.DE2_D2, first);
    dump_field(out, "DE2_MS", bl.DE2_MS, first);

    dump_field(out, "CFM", bl.CFM, first);
    dump_field(out, "CFM_HKA", bl.CFM_HKA, first);
    dump_field(out, "CFM_MA", bl.CFM_MA, first);
    dump_field(out, "CFM_MS", bl.CFM_MS, first);
    dump_field(out, "CFM_RE", bl.CFM_RE, first);
    dump_field(out, "CFM_RTA", bl.CFM_RTA, first);
    dump_field(out, "CFM_U1", bl.CFM_U1, first);
    dump_field(out, "CFM_T1", bl.CFM_T1, first);
    dump_field(out, "CFM_D1", bl.CFM_D1, first);
    dump_field(out, "CFM_U2", bl.CFM_U2, first);
    dump_field(out, "CFM_T2", bl.CFM_T2, first);
    dump_field(out, "CFM_D2", bl.CFM_D2, first);
    dump_field(out, "XT", bl.XT, first);
    dump_field(out, "XT_A1", bl.XT_A1, first);
    dump_field(out, "XT_A2", bl.XT_A2, first);
    dump_field(out, "XT_MS", bl.XT_MS, first);
    dump_field(out, "XT_RE", bl.XT_RE, first);
    dump_field(out, "XT_XF", bl.XT_XF, first);
    dump_field(out, "XT_X1", bl.XT_X1, first);
    dump_field(out, "XT_T1", bl.XT_T1, first);
    dump_field(out, "XT_D1", bl.XT_D1, first);
    dump_field(out, "XT_U1", bl.XT_U1, first);
    dump_field(out, "XT_X2", bl.XT_X2, first);
    dump_field(out, "XT_T2", bl.XT_T2, first);
    dump_field(out, "XT_D2", bl.XT_D2, first);
    dump_field(out, "XT_U2", bl.XT_U2, first);
    dump_field(out, "DWTE", bl.DWTE, first);
    dump_field(out, "QINFBL", bl.QINFBL, first);
    dump_field(out, "TKBL", bl.TKBL, first);
    dump_field(out, "TKBL_MS", bl.TKBL_MS, first);
    dump_field(out, "RSTBL", bl.RSTBL, first);
    dump_field(out, "RSTBL_MS", bl.RSTBL_MS, first);
    dump_field(out, "HSTINV", bl.HSTINV, first);
    dump_field(out, "HSTINV_MS", bl.HSTINV_MS, first);
    dump_field(out, "REYBL", bl.REYBL, first);
    dump_field(out, "REYBL_MS", bl.REYBL_MS, first);
    dump_field(out, "REYBL_RE", bl.REYBL_RE, first);
    dump_field(out, "GAMBL", bl.GAMBL, first);
    dump_field(out, "GM1BL", bl.GM1BL, first);
    dump_field(out, "HVRAT", bl.HVRAT, first);
    dump_field(out, "BULE", bl.BULE, first);
    dump_field(out, "XIFORC", bl.XIFORC, first);
    dump_field(out, "AMCRIT", bl.AMCRIT, first);
    out << '}';
}

}  // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: compare_xbl <state_json> [--out path] [--mode setbl|mrchue|mrchdu|viscal] [--niter n]\n";
        return 2;
    }

    std::string state_path;
    std::string out_path;
    std::string mode = "setbl";
    int niter = 1;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--out" && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            mode = argv[++i];
        } else if (arg == "--niter" && i + 1 < argc) {
            niter = std::stoi(argv[++i]);
        } else if (state_path.empty()) {
            state_path = arg;
        }
    }

    if (state_path.empty()) {
        std::cerr << "Missing state_json path\n";
        return 2;
    }

    std::ifstream handle(state_path);
    if (!handle) {
        std::cerr << "Failed to open " << state_path << "\n";
        return 2;
    }

    std::stringstream buffer;
    buffer << handle.rdbuf();
    Json data = parse_json(buffer.str());
    if (data.type != Json::Type::kObject) {
        throw std::runtime_error("Root is not an object");
    }

    XFoilState ctx;
    XBlState bl;

    std::unordered_map<std::string, std::function<void(const Json &)>> ctx_setters;
    std::unordered_map<std::string, std::function<void(const Json &)>> bl_setters;
    populate_ctx_setters(ctx, ctx_setters);
    populate_bl_setters(bl, bl_setters);

    const auto ctx_it = data.object.find("ctx");
    if (ctx_it != data.object.end()) {
        const Json &ctx_json = ctx_it->second;
        if (ctx_json.type == Json::Type::kObject) {
            for (const auto &kv : ctx_json.object) {
                const auto setter = ctx_setters.find(kv.first);
                if (setter != ctx_setters.end()) {
                    setter->second(kv.second);
                }
            }
        }
    }

    const auto bl_it = data.object.find("bl");
    if (bl_it != data.object.end()) {
        const Json &bl_json = bl_it->second;
        if (bl_json.type == Json::Type::kObject) {
            for (const auto &kv : bl_json.object) {
                const auto setter = bl_setters.find(kv.first);
                if (setter != bl_setters.end()) {
                    setter->second(kv.second);
                }
            }
        }
    }

    if (mode == "setbl") {
        setbl(ctx, bl);
    } else if (mode == "mrchue") {
        mrchue(ctx, bl);
    } else if (mode == "mrchdu") {
        mrchdu(ctx, bl);
    } else if (mode == "viscal") {
        viscal(ctx, bl, niter);
    } else {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 2;
    }

    std::ostringstream out_buf;
    out_buf << '{';
    out_buf << "\"ctx\":";
    dump_ctx(out_buf, ctx);
    out_buf << ",\"bl\":";
    dump_bl(out_buf, bl);
    out_buf << '}';

    if (!out_path.empty()) {
        std::ofstream out_file(out_path);
        out_file << out_buf.str() << "\n";
    } else {
        std::cout << out_buf.str() << "\n";
    }

    return 0;
}
