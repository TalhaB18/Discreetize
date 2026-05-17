#include "supabase_client.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <sstream>
#include <fstream>

namespace cfd {

std::string SupabaseClient::url_;
std::string SupabaseClient::service_key_;
std::string SupabaseClient::anon_key_;
bool        SupabaseClient::configured_ = false;

bool SupabaseClient::init() {
    const char* u = std::getenv("SUPABASE_URL");
    const char* sk = std::getenv("SUPABASE_SERVICE_KEY");
    const char* ak = std::getenv("SUPABASE_ANON_KEY");
    if(!u || !*u) { configured_ = false; return false; }
    url_ = u;
    if(sk) service_key_ = sk;
    if(ak) anon_key_ = ak;
    // Remove trailing slash
    while(!url_.empty() && url_.back()=='/') url_.pop_back();
    configured_ = true;
    return true;
}

std::string SupabaseClient::getPublicUrl(const std::string& bucket,
                                          const std::string& storage_path) {
    if(!configured_) return "";
    return url_ + "/storage/v1/object/public/" + bucket + "/" + storage_path;
}

std::string SupabaseClient::uploadFile(const std::string& bucket,
                                        const std::string& storage_path,
                                        const std::string& data,
                                        const std::string& content_type) {
    if(!configured_ || service_key_.empty()) return "";

    // Write data to temp file
    std::string tmp = "/tmp/sb_upload_" + std::to_string(
        static_cast<unsigned long>(std::time(nullptr)) ^ static_cast<unsigned long>(rand())) + ".bin";
    {
        std::ofstream f(tmp, std::ios::binary);
        if(!f) return "";
        f.write(data.data(), static_cast<std::streamsize>(data.size()));
    }

    std::string endpoint = url_ + "/storage/v1/object/" + bucket + "/" + storage_path + "?upsert=true";
    std::string auth = "Authorization: Bearer " + service_key_;
    std::string body_arg = "@" + tmp;

    std::ostringstream cmd;
    cmd << "curl -s -o /dev/null -w \"%{http_code}\" -X PUT"
        << " -H \"" << auth << "\""
        << " -H \"Content-Type: " << content_type << "\""
        << " --data-binary '" << body_arg << "'"
        << " \"" << endpoint << "\"";

    std::string code = httpCurl("PUT_FILE", tmp, data, content_type, auth);
    std::remove(tmp.c_str());

    if(code.empty()) return "";
    return getPublicUrl(bucket, storage_path);
}

std::string SupabaseClient::insertRow(const std::string& table,
                                       const std::string& json_body) {
    if(!configured_ || service_key_.empty()) return "";
    std::string endpoint = url_ + "/rest/v1/" + table;
    std::string auth = "Authorization: Bearer " + service_key_;
    std::string resp = httpCurl("POST", endpoint, json_body, "application/json", auth);
    // Extract id from response JSON
    const std::string key = "\"id\"";
    size_t kp = resp.find(key);
    if(kp == std::string::npos) return "";
    size_t cp = resp.find(':', kp);
    if(cp == std::string::npos) return "";
    size_t qs = resp.find('"', cp+1);
    if(qs == std::string::npos) return "";
    size_t qe = resp.find('"', qs+1);
    if(qe == std::string::npos) return "";
    return resp.substr(qs+1, qe-qs-1);
}

bool SupabaseClient::updateRow(const std::string& table, const std::string& id,
                                const std::string& json_body) {
    if(!configured_ || service_key_.empty()) return false;
    std::string endpoint = url_ + "/rest/v1/" + table + "?id=eq." + id;
    std::string auth = "Authorization: Bearer " + service_key_;
    std::string resp = httpCurl("PATCH", endpoint, json_body, "application/json", auth);
    return !resp.empty();
}

std::string SupabaseClient::downloadFile(const std::string& bucket,
                                          const std::string& storage_path) {
    if(!configured_ || service_key_.empty()) return "";

    std::string endpoint = url_ + "/storage/v1/object/" + bucket + "/" + storage_path;
    std::string auth = "Authorization: Bearer " + service_key_;

    // Write output to a temp file to handle binary data safely
    std::string tmp = "/tmp/sb_dl_" + std::to_string(static_cast<unsigned long>(std::time(nullptr))
                       ^ static_cast<unsigned long>(rand())) + ".bin";

    std::ostringstream cmd;
    cmd << "curl -s -o \"" << tmp << "\""
        << " -H \"" << auth << "\""
        << " \"" << endpoint << "\" 2>/dev/null";

    FILE* pipe = popen(cmd.str().c_str(), "r");
    if(pipe) pclose(pipe);

    // Read the downloaded file
    std::ifstream f(tmp, std::ios::binary);
    if(!f.is_open()) { std::remove(tmp.c_str()); return ""; }
    std::ostringstream ss;
    ss << f.rdbuf();
    std::remove(tmp.c_str());
    return ss.str();
}

std::string SupabaseClient::queryRows(const std::string& table,
                                       const std::string& filter) {
    if(!configured_ || service_key_.empty()) return "";

    std::string endpoint = url_ + "/rest/v1/" + table;
    if(!filter.empty()) endpoint += "?" + filter;

    std::string auth = "Authorization: Bearer " + service_key_;

    // Build curl command (GET — pass empty body)
    std::string tmp = "/tmp/sb_qbody_" + std::to_string(rand()) + ".json";
    {
        std::ofstream f(tmp);
        f << "";
    }

    std::ostringstream cmd;
    cmd << "curl -s -X GET"
        << " -H \"" << auth << "\""
        << " -H \"Content-Type: application/json\""
        << " -H \"Prefer: return=representation\""
        << " \"" << endpoint << "\" 2>/dev/null";

    std::string result;
    FILE* pipe = popen(cmd.str().c_str(), "r");
    if(pipe) {
        char buf[4096];
        while(fgets(buf, sizeof(buf), pipe)) result += buf;
        pclose(pipe);
    }
    std::remove(tmp.c_str());
    return result;
}

bool SupabaseClient::deleteRow(const std::string& table, const std::string& filter) {
    if(!configured_ || service_key_.empty()) return false;

    std::string endpoint = url_ + "/rest/v1/" + table;
    if(!filter.empty()) endpoint += "?" + filter;

    std::string auth = "Authorization: Bearer " + service_key_;

    std::ostringstream cmd;
    cmd << "curl -s -o /dev/null -w \"%{http_code}\" -X DELETE"
        << " -H \"" << auth << "\""
        << " -H \"Content-Type: application/json\""
        << " -H \"Prefer: return=minimal\""
        << " \"" << endpoint << "\" 2>/dev/null";

    std::string code;
    FILE* pipe = popen(cmd.str().c_str(), "r");
    if(pipe) {
        char buf[64];
        while(fgets(buf, sizeof(buf), pipe)) code += buf;
        pclose(pipe);
    }
    // 2xx = success
    if(code.size() >= 3) {
        int http_code = 0;
        try { http_code = std::stoi(code.substr(0,3)); } catch(...) {}
        return http_code >= 200 && http_code < 300;
    }
    return false;
}

std::string SupabaseClient::httpCurl(const std::string& method,
                                      const std::string& endpoint,
                                      const std::string& body,
                                      const std::string& content_type,
                                      const std::string& auth_header) {
    // Write body to temp file
    std::string tmp = "/tmp/sb_body_" + std::to_string(rand()) + ".json";
    {
        std::ofstream f(tmp);
        f << body;
    }

    std::ostringstream cmd;
    cmd << "curl -s -X " << method
        << " -H \"" << auth_header << "\""
        << " -H \"Content-Type: " << content_type << "\""
        << " -H \"Prefer: return=representation\""
        << " --data-binary @" << tmp
        << " \"" << endpoint << "\" 2>/dev/null";

    std::string result;
    FILE* pipe = popen(cmd.str().c_str(), "r");
    if(pipe){
        char buf[4096];
        while(fgets(buf, sizeof(buf), pipe)) result += buf;
        pclose(pipe);
    }
    std::remove(tmp.c_str());
    return result;
}

} // namespace cfd
