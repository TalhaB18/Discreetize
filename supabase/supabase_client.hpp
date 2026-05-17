#pragma once
#include <string>

namespace cfd {

// Supabase integration. Configure via env vars:
//   SUPABASE_URL         https://<project>.supabase.co
//   SUPABASE_SERVICE_KEY server-side service role key
//   SUPABASE_ANON_KEY    public anon key
//
// All HTTP calls use system curl (no libcurl dependency).
// When not configured, all methods return "" gracefully.
class SupabaseClient {
public:
    static bool init();
    static bool isConfigured() { return configured_; }
    static const std::string& url() { return url_; }

    // Upload bytes to Supabase Storage. Returns public URL or "".
    static std::string uploadFile(const std::string& bucket,
                                  const std::string& storage_path,
                                  const std::string& data,
                                  const std::string& content_type = "application/octet-stream");

    // Insert JSON row into a PostgREST table. Returns created id or "".
    static std::string insertRow(const std::string& table, const std::string& json_body);

    // PATCH update a row. Returns true on success.
    static bool updateRow(const std::string& table, const std::string& id,
                          const std::string& json_body);

    // Construct public URL without HTTP request.
    static std::string getPublicUrl(const std::string& bucket,
                                    const std::string& storage_path);

    // Download a file from Supabase Storage. Returns file bytes or "" on error.
    static std::string downloadFile(const std::string& bucket,
                                    const std::string& storage_path);

    // Query rows from a table with an optional filter (e.g. "id=eq.abc123").
    // Returns raw JSON response string.
    static std::string queryRows(const std::string& table,
                                 const std::string& filter = "");

    // Delete a row (or rows) matching filter. Returns true on success.
    static bool deleteRow(const std::string& table, const std::string& filter);

private:
    static std::string url_;
    static std::string service_key_;
    static std::string anon_key_;
    static bool        configured_;

    static std::string httpCurl(const std::string& method,
                                const std::string& endpoint,
                                const std::string& body,
                                const std::string& content_type,
                                const std::string& auth_header);
};

} // namespace cfd
